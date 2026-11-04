"""Kaggle - Challenge 5 phase 4. Score every design against all three states.

WHAT THIS IS FOR. The competition ranks on two ratios before it looks at strength:
KD(wild type)/KD(G12D) and KD(peptide-free)/KD(G12D). So the quantity we need to
predict is a MARGIN between near-identical complexes, not an affinity.

    g12d   HLA + b2m + VVVGADGVGK     must bind
    wt     HLA + b2m + VVVGAGGVGK     must not bind   (one residue different)
    apo    HLA + b2m                  must not bind   (no peptide)

THE CONTROL THAT MAKES THE MARGIN MEANINGFUL. AlphaFold's confidence varies between its
own five parameter sets for identical input. If i_ptm differs by 0.04 between model 1
and model 2 on the SAME state, then a 0.03 difference between g12d and wt is noise with
a story attached. So every state is predicted with several models and the spread is
reported beside the margin. A margin smaller than the spread is reported as not
measurable, not as a small effect.

This is the control that four previous challenges did not have, and it is the difference
between a ranking and a ranking we can defend.

WHAT THIS CANNOT DO, STATED HERE SO IT REACHES THE METHODS TEXT
  - The wild-type state carries the mutant's backbone with Asp6 cut back to Gly. That
    isolates the substitution, at the cost of not being a real wild-type structure.
  - The apo state is a holo groove with the peptide deleted. Real peptide-free HLA is
    unstable and partly collapsed, so this is an easier test than the assay performs.
  - AlphaFold is not known to predict specificity margins between near-identical
    complexes. A predicted margin is a hypothesis, not a measurement.

PAE is saved for every prediction, so ipSAE can be computed afterwards on a processor.
"""
import glob
import json
import os
import subprocess
import sys
import time

W = "/kaggle/working"
MODELS = (0, 1)                 # two parameter sets per state: the noise floor
CONTACT = 4.5
MUT_POS = 6
PEP_G12D = "VVVGADGVGK"
PEP_WT = "VVVGAGGVGK"
BACKBONE = ("N", "CA", "C", "O", "OXT")
BUDGET_MIN = 240


def sh(c):
    subprocess.run(c, shell=True, check=True)


try:
    import colabdesign  # noqa: F401
except ImportError:
    sh("pip -q install git+https://github.com/sokrypton/ColabDesign.git@v1.1.1 biopython")

import copy  # noqa: E402

import numpy as np  # noqa: E402
from Bio import PDB  # noqa: E402

if not os.path.exists("%s/params/params_model_1_ptm.npz" % W):
    os.makedirs("%s/params" % W, exist_ok=True)
    sh("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
       "alphafold_params_2022-12-06.tar" % W)
    sh("tar -xf %s/af2.tar -C %s/params" % (W, W))
    os.remove("%s/af2.tar" % W)
if len(glob.glob("%s/params/params_model_*.npz" % W)) == 0:
    raise SystemExit("no AlphaFold parameters - refusing to start")

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def residues(ch):
    out, seen = [], set()
    for r in ch:
        if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def seq(rs):
    return "".join(T3[r.get_resname()] for r in rs)


def heavy(r):
    return np.array([a.coord for a in r if a.element != "H"])


class Sel(PDB.Select):
    def __init__(self, keep):
        self.keep = {c: set(n) for c, n in keep.items()}

    def accept_model(self, m):
        return m.id == 0

    def accept_chain(self, c):
        return c.id in self.keep

    def accept_residue(self, r):
        return r.id[0] == " " and r.id[1] in self.keep[r.get_parent().id]

    def accept_atom(self, a):
        return (not a.is_disordered()) or a.get_altloc() in (" ", "A")


def build_states():
    """Three structures sharing one background, differing only as intended."""
    import urllib.request
    for pid in ("7OW6", "7OW5"):
        f = "%s/%s.cif" % (W, pid)
        if not os.path.exists(f):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/%s.cif" % pid, f)
    p = PDB.MMCIFParser(QUIET=True)
    m6 = p.get_structure("6", "%s/7OW6.cif" % W)[0]
    m5 = p.get_structure("5", "%s/7OW5.cif" % W)[0]
    hv_all, b2m_all, pep = residues(m6["A"]), residues(m6["B"]), residues(m6["C"])
    if seq(pep) != PEP_G12D:
        raise SystemExit("7OW6 chain C is %s" % seq(pep))
    hv_ok = {r.id[1] for r in hv_all} & {r.id[1] for r in residues(m5["A"])}
    b2m_ok = {r.id[1] for r in b2m_all} & {r.id[1] for r in residues(m5["B"])}
    holo = {"A": sorted(hv_ok), "B": sorted(b2m_ok), "C": [r.id[1] for r in pep]}
    apo = {k: v for k, v in holo.items() if k != "C"}
    mut_resnum = pep[MUT_POS - 1].id[1]

    io = PDB.PDBIO()
    io.set_structure(m6)
    io.save("%s/state_g12d.pdb" % W, Sel(holo))
    io.save("%s/state_apo.pdb" % W, Sel(apo))

    mut = copy.deepcopy(m6)
    tgt = mut["C"][mut_resnum]
    for nm in [a.get_name() for a in tgt if a.get_name() not in BACKBONE]:
        tgt.detach_child(nm)
    tgt.resname = "GLY"
    io.set_structure(mut)
    io.save("%s/state_wt.pdb" % W, Sel(holo))

    states = {}
    for tag, chains in (("g12d", "A,B,C"), ("wt", "A,B,C"), ("apo", "A,B")):
        path = "%s/state_%s.pdb" % (W, tag)
        m = PDB.PDBParser(QUIET=True).get_structure("s", path)[0]
        ch = {c.id: residues(c) for c in m}
        pseq = seq(ch["C"]) if "C" in ch else None
        states[tag] = dict(path=path, chains=chains, peptide=pseq,
                           n=sum(len(v) for v in ch.values()))
    if states["g12d"]["peptide"] != PEP_G12D:
        raise SystemExit("g12d peptide is %s" % states["g12d"]["peptide"])
    if states["wt"]["peptide"] != PEP_WT:
        raise SystemExit("wt peptide is %s - the mutation was not applied"
                         % states["wt"]["peptide"])
    if states["apo"]["peptide"] is not None:
        raise SystemExit("apo state still carries a peptide")
    if states["g12d"]["n"] != states["wt"]["n"]:
        raise SystemExit("g12d and wt are not the same size: %d vs %d"
                         % (states["g12d"]["n"], states["wt"]["n"]))
    for tag, s in states.items():
        print("state %-5s %3d residues  peptide %s" % (tag, s["n"], s["peptide"]))
    return states


def contacts(path, state_tag, pep_seq):
    """Binder contacts on the peptide, the mutation, and the groove."""
    m = PDB.PDBParser(QUIET=True).get_structure("d", path)[0]
    ch = {c.id: residues(c) for c in m}
    sizes = {c: len(rs) for c, rs in ch.items()}
    # the target chains are merged by save_pdb; the binder is whatever is left
    if len(ch) != 2:
        raise RuntimeError("expected two chains after save_pdb, got %s" % sizes)
    tgt = max(ch, key=lambda c: len(ch[c]))
    binder = [c for c in ch if c != tgt][0]
    b = np.concatenate([heavy(r) for r in ch[binder]])

    def touch(rs):
        return sum(1 for r in rs
                   if np.linalg.norm(heavy(r)[:, None, :] - b[None, :, :],
                                     axis=-1).min() <= CONTACT)

    if pep_seq is None:
        return dict(on_mutation=0, on_peptide=0, on_groove=touch(ch[tgt]))
    rs = ch[tgt]
    if seq(rs)[-len(pep_seq):] != pep_seq:
        raise RuntimeError("%s: target chain does not end in the peptide %s"
                           % (state_tag, pep_seq))
    pep_rs, groove_rs = rs[-len(pep_seq):], rs[:-len(pep_seq)]
    return dict(on_mutation=touch([pep_rs[MUT_POS - 1]]), on_peptide=touch(pep_rs),
                on_groove=touch(groove_rs))


def main():
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
    states = build_states()

    src = None
    for cand in ("%s/pmhc_gen.json" % W,
                 *sorted(glob.glob("/kaggle/input/**/pmhc_gen.json", recursive=True))):
        if os.path.exists(cand):
            src = cand
            break
    if src is None:
        raise SystemExit("no pmhc_gen.json - attach the generation run's output")
    designs = json.load(open(src))
    designs = [d for d in designs if d.get("sequence")]
    print("\n%d designs from %s\n" % (len(designs), src))
    if not designs:
        raise SystemExit("generation produced no sequences - nothing to screen")

    from colabdesign import clear_mem, mk_afdesign_model

    t0 = time.time()
    out_path = "%s/pmhc_screen.jsonl" % W
    done = 0
    for d in designs:
        if (time.time() - t0) / 60.0 >= BUDGET_MIN:
            print("budget reached, stopping after %d designs" % done)
            break
        rec = dict(name=d["name"], arm=d.get("arm"), sequence=d["sequence"],
                   length=len(d["sequence"]), states={})
        try:
            for tag, st in states.items():
                per_model = []
                for mi in MODELS:
                    clear_mem()
                    model = mk_afdesign_model(protocol="binder", data_dir=W)
                    model.prep_inputs(pdb_filename=st["path"], chain=st["chains"],
                                      binder_len=len(d["sequence"]))
                    model.predict(seq=d["sequence"], models=[mi], verbose=False)
                    lg = model.aux.get("log", {})
                    pdb_out = "%s/%s_%s_m%d.pdb" % (W, d["name"], tag, mi)
                    model.save_pdb(pdb_out)
                    pae = model.aux.get("pae")
                    if pae is not None:
                        np.save("%s/%s_%s_m%d_pae.npy" % (W, d["name"], tag, mi),
                                np.asarray(pae))
                    c = contacts(pdb_out, tag, st["peptide"])
                    per_model.append(dict(model=mi,
                                          i_ptm=round(float(lg.get("i_ptm", 0)), 4),
                                          ptm=round(float(lg.get("ptm", 0)), 4),
                                          plddt=round(float(lg.get("plddt", 0)), 4),
                                          pae_saved=pae is not None, **c))
                rec["states"][tag] = dict(
                    per_model=per_model,
                    i_ptm=round(float(np.mean([x["i_ptm"] for x in per_model])), 4),
                    i_ptm_spread=round(float(np.ptp([x["i_ptm"] for x in per_model])), 4),
                    on_mutation=max(x["on_mutation"] for x in per_model),
                    on_peptide=max(x["on_peptide"] for x in per_model),
                    on_groove=max(x["on_groove"] for x in per_model))
        except Exception as e:
            print("%s FAILED %s %s" % (d["name"], type(e).__name__, str(e)[:160]))
            continue

        g, w, a = rec["states"]["g12d"], rec["states"]["wt"], rec["states"]["apo"]
        noise = max(g["i_ptm_spread"], w["i_ptm_spread"], a["i_ptm_spread"])
        rec["margin_wt"] = round(g["i_ptm"] - w["i_ptm"], 4)
        rec["margin_apo"] = round(g["i_ptm"] - a["i_ptm"], 4)
        rec["noise"] = noise
        rec["measurable"] = bool(min(rec["margin_wt"], rec["margin_apo"]) > noise)
        with open(out_path, "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        done += 1
        print("%-22s g12d %.3f wt %.3f apo %.3f | margin wt %+.3f apo %+.3f | "
              "noise %.3f | %s | on-mut g12d %d [%.1f min]"
              % (d["name"], g["i_ptm"], w["i_ptm"], a["i_ptm"], rec["margin_wt"],
                 rec["margin_apo"], noise,
                 "MEASURABLE" if rec["measurable"] else "within noise",
                 g["on_mutation"], (time.time() - t0) / 60.0))

    records = [json.loads(l) for l in open(out_path)] if os.path.exists(out_path) else []
    json.dump(records, open("%s/pmhc_screen.json" % W, "w"), indent=0)

    print("\n==================== COUNTER-SCREEN RESULT ====================")
    if not records:
        raise SystemExit("nothing was screened")
    noise = float(np.median([r["noise"] for r in records]))
    print("Median model-to-model spread in i_ptm on identical input: %.3f" % noise)
    print("Any margin below that is noise with a story attached.\n")
    good = [r for r in records if r["measurable"]]
    print("%d of %d designs have BOTH margins above their own noise floor."
          % (len(good), len(records)))
    print("\n%-22s %7s %7s %7s %8s %8s %6s"
          % ("design", "g12d", "wt", "apo", "marg_wt", "marg_apo", "on-mut"))
    for r in sorted(records, key=lambda r: -min(r["margin_wt"], r["margin_apo"]))[:25]:
        s = r["states"]
        print("%-22s %7.3f %7.3f %7.3f %+8.3f %+8.3f %6d%s"
              % (r["name"], s["g12d"]["i_ptm"], s["wt"]["i_ptm"], s["apo"]["i_ptm"],
                 r["margin_wt"], r["margin_apo"], s["g12d"]["on_mutation"],
                 "  *" if r["measurable"] else ""))
    print("\n===== COPY FROM HERE =====")
    for r in records:
        print(json.dumps(r))
    print("===== COPY TO HERE =====")


main()
