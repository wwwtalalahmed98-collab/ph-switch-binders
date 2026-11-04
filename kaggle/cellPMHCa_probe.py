"""Kaggle - Challenge 5 probe. Does the tight hotspot set actually reach the mutation?

WHAT THIS IS FOR. Phase 1 chose peptide positions 5,6,7 over the obvious set of 4-9, on
the argument that Gly7 at 73% exposure is six times easier to reach than Asp6 at 11%, so
a term spread over six positions would be satisfied by Gly7 alone and never read the
mutation. That is an argument, not a measurement. This run measures it, on two arms, for
about forty minutes, before any long run is committed.

    arm A   hotspots C5,C6,C7       num_pos 3    the chosen set
    arm B   hotspots C4..C9         num_pos 6    the obvious set

If arm A does not beat arm B on contacts with residue 6 specifically, the phase 1
reasoning is wrong and the long run must not use it.

THE METRIC. on_mutation - binder heavy atoms within 4.5 A of PEPTIDE RESIDUE 6 - and
nothing else, because on TfR1 a round was reported as a success on a metric that counted
contacts anywhere near the target. A design with twenty peptide contacts and none on the
mutation reads nothing the wild type lacks.

THREE THINGS THAT WENT WRONG BEFORE, GUARDED HERE
  1. i_con's default cutoff is 21.69 A, which passes 63 of 64 distance bins and accepts
     "contacts" out to 21.4 A. TfR1 round 1 aimed at nothing because of it. This sets
     8.0 and REFUSES TO RUN if the setting did not take.
  2. save_pdb renumbers residues from 1. Two TfR1 rounds were measured against residue
     numbers that no longer existed and silently returned nothing. The target chain here
     is found by reproducing its sequence residue for residue.
  3. We have never once computed ipSAE, in four challenges, for one reason: it needs the
     PAE matrix and we never saved it. This saves it per design, from the first run.

THE TARGET IS REBUILT, NOT UPLOADED. The notebook fetches 7OW6 and 7OW5 and rebuilds the
20 A crop by the same rule as pmhc/build_states.py, then checks the result against the
sequence that rule produced locally. If a single residue differs it stops. That is
cheaper than carrying 41 kB of base64 and it verifies itself.
"""
import json
import os
import time

import numpy as np

W = "/kaggle/working"
BUDGET_MIN = 40
LEN_MIN, LEN_MAX = 65, 90
SEED0 = 7100

MUT_POS = 6
CONTACT = 4.5
PEP_G12D = "VVVGADGVGK"
DESIGN_RADIUS = 20.0
MIN_SEGMENT = 3

EXPECTED_A = ("SHSMRYFYTSVSRPGGEPRFIAVGYVDDTQFVRFDSDAASQRMEPRAPWIEQEGPEYWDQETRNVKAQSQ"
              "TDRVDLGTLRGYYNQSEDGSHTIQIMYGCDVGPDGRFLRGYRQDAYDGKDYIALNEDLRSWTAADMAAQI"
              "TKRKWEAAHAAEQQRAYLEGRCVEWLRRYLENGKETLQR")

ARMS = {
    "A_tight": dict(hotspot=[5, 6, 7], num_pos=3),
    "B_broad": dict(hotspot=[4, 5, 6, 7, 8, 9], num_pos=6),
}

T3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
      "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
      "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
      "TYR": "Y", "VAL": "V"}


def residues(chain):
    """Ordered, de-duplicated standard residues, as a writer emits them."""
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T3 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def seq(rs):
    return "".join(T3[r.get_resname()] for r in rs)


def heavy(r):
    return np.array([a.coord for a in r if a.element != "H"])


def build_target(path):
    """Rebuild the 20 A crop and prove it is the one phase 0 produced."""
    import urllib.request
    from Bio import PDB

    for pid in ("7OW6", "7OW5"):
        f = "%s/%s.cif" % (W, pid)
        if not os.path.exists(f):
            urllib.request.urlretrieve(
                "https://files.rcsb.org/download/%s.cif" % pid, f)
    p = PDB.MMCIFParser(QUIET=True)
    m6 = p.get_structure("6", "%s/7OW6.cif" % W)[0]
    m5 = p.get_structure("5", "%s/7OW5.cif" % W)[0]

    hv, pep = residues(m6["A"]), residues(m6["C"])
    if seq(pep) != PEP_G12D:
        raise SystemExit("7OW6 chain C is %s, expected %s" % (seq(pep), PEP_G12D))
    allowed = {r.id[1] for r in hv} & {r.id[1] for r in residues(m5["A"])}

    pep_atoms = np.concatenate([heavy(r) for r in pep])
    near = [r.id[1] for r in hv if r.id[1] in allowed
            and np.linalg.norm(heavy(r)[:, None, :] - pep_atoms[None, :, :],
                               axis=-1).min() <= DESIGN_RADIUS]
    runs, cur = [], []
    for n in sorted(near):
        if cur and n == cur[-1] + 1:
            cur.append(n)
        else:
            if len(cur) >= MIN_SEGMENT:
                runs += cur
            cur = [n]
    if len(cur) >= MIN_SEGMENT:
        runs += cur

    keep = {"A": set(runs), "C": {r.id[1] for r in pep}}

    class Sel(PDB.Select):
        def accept_model(self, m):
            return m.id == 0

        def accept_chain(self, c):
            return c.id in keep

        def accept_residue(self, r):
            return r.id[0] == " " and r.id[1] in keep[r.get_parent().id]

        def accept_atom(self, a):
            return (not a.is_disordered()) or a.get_altloc() in (" ", "A")

    io = PDB.PDBIO()
    io.set_structure(m6)
    io.save(path, Sel())

    chk = PDB.PDBParser(QUIET=True).get_structure("t", path)[0]
    got_a, got_c = seq(residues(chk["A"])), seq(residues(chk["C"]))
    if got_a != EXPECTED_A:
        raise SystemExit("rebuilt chain A does not match phase 0:\n  got %d aa\n  want "
                         "%d aa" % (len(got_a), len(EXPECTED_A)))
    if got_c != PEP_G12D:
        raise SystemExit("rebuilt peptide is %s" % got_c)
    print("target rebuilt and verified: chain A %d aa, peptide %s" % (len(got_a), got_c))
    return path


def measure(pdb_path, crop_a, crop_c):
    """Contacts, counted against the peptide by sequence match - never by numbering.

    save_pdb MERGES the two target chains into one. Version 2 of this notebook died
    here with chains {'A': 189, 'B': 77}: 179 groove residues and the 10-residue
    peptide concatenated, plus the binder. So the target is found as one chain whose
    sequence is the two crops joined, and the peptide is located by WHERE the join
    falls - checked both ways round, because which half comes first is save_pdb's
    business, not something to assume.
    """
    from Bio import PDB
    m = PDB.PDBParser(QUIET=True).get_structure("d", pdb_path)[0]
    chains = {c.id: residues(c) for c in m}

    tgt = pep_rs = None
    for cid, rs in chains.items():
        s = seq(rs)
        if s == crop_a + crop_c:
            tgt, pep_rs = cid, rs[len(crop_a):]
            break
        if s == crop_c + crop_a:
            tgt, pep_rs = cid, rs[:len(crop_c)]
            break
    if tgt is None:
        # fall back to the separate-chain layout, in case save_pdb ever changes
        a = next((c for c, rs in chains.items() if seq(rs) == crop_a), None)
        c_ = next((c for c, rs in chains.items()
                   if c != a and seq(rs) == crop_c), None)
        if a is None or c_ is None:
            raise SystemExit("no chain reproduces the crop; chains are %s"
                             % {c: len(rs) for c, rs in chains.items()})
        tgt, pep_rs = a, chains[c_]
        groove_rs = chains[a]
        binder = [c for c in chains if c not in (a, c_)]
    else:
        groove_rs = [r for r in chains[tgt] if r not in pep_rs]
        binder = [c for c in chains if c != tgt]

    if seq(pep_rs) != crop_c:
        raise SystemExit("located peptide reads %s, expected %s"
                         % (seq(pep_rs), crop_c))
    if len(binder) != 1:
        raise SystemExit("expected one binder chain, got %s" % binder)
    b_atoms = np.concatenate([heavy(r) for r in chains[binder[0]]])

    def touch(rs):
        n = 0
        for r in rs:
            if np.linalg.norm(heavy(r)[:, None, :] - b_atoms[None, :, :],
                              axis=-1).min() <= CONTACT:
                n += 1
        return n

    on_mut = touch([pep_rs[MUT_POS - 1]])
    return dict(on_mutation=on_mut, on_peptide=touch(pep_rs),
                on_groove=touch(groove_rs),
                mutation_contacts=int(np.sum(
                    np.linalg.norm(heavy(pep_rs[MUT_POS - 1])[:, None, :]
                                   - b_atoms[None, :, :], axis=-1) <= CONTACT)))


def fetch_params():
    """AlphaFold weights. Without these ColabDesign builds a model with no parameters
    and dies at the first optimisation step with 'no model params defined' - which is
    exactly how version 1 of this notebook failed, 47 seconds in."""
    import glob
    import subprocess

    pdir = "%s/params" % W
    if not os.path.exists("%s/params_model_1_ptm.npz" % pdir):
        os.makedirs(pdir, exist_ok=True)
        subprocess.run("wget -q -O %s/af2.tar https://storage.googleapis.com/alphafold/"
                       "alphafold_params_2022-12-06.tar" % W, shell=True, check=True)
        subprocess.run("tar -xf %s/af2.tar -C %s" % (W, pdir), shell=True, check=True)
        os.remove("%s/af2.tar" % W)
    n = len(glob.glob("%s/params_model_*.npz" % pdir))
    print("params: %d files" % n)
    if n == 0:
        raise SystemExit("no AlphaFold parameters - refusing to start, the run would "
                         "fail at the first design step and waste the budget")


def main():
    os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
    os.environ["XLA_PYTHON_CLIENT_ALLOCATOR"] = "platform"
    fetch_params()
    target = build_target("%s/target_design.pdb" % W)

    from colabdesign import clear_mem, mk_afdesign_model

    crop_a, crop_c = EXPECTED_A, PEP_G12D
    t0 = time.time()
    records = []
    n = 0
    while (time.time() - t0) / 60.0 < BUDGET_MIN:
        for arm, cfg in ARMS.items():
            if (time.time() - t0) / 60.0 >= BUDGET_MIN:
                break
            seed = SEED0 + n
            rng = np.random.RandomState(seed)
            L = int(rng.randint(LEN_MIN, LEN_MAX + 1))
            hot = ",".join("C%d" % p for p in cfg["hotspot"])
            clear_mem()
            model = mk_afdesign_model(protocol="binder", data_dir=W)
            model.prep_inputs(pdb_filename=target, chain="A,C", binder_len=L,
                              hotspot=hot)

            before = dict(model.opt["i_con"])
            if "cutoff" not in before:
                raise SystemExit("this colabdesign has no i_con cutoff: %s" % before)
            model.set_opt(i_con=dict(cutoff=8.0, num=1, num_pos=cfg["num_pos"]),
                          set_defaults=True)
            now = dict(model.opt["i_con"])
            if abs(float(now["cutoff"]) - 8.0) > 0.01:
                raise SystemExit("i_con cutoff is still %s - refusing to repeat TfR1 "
                                 "round 1" % now["cutoff"])
            if n == 0:
                print("i_con %s -> %s\n" % (before, now))

            model.restart(seed=seed)
            model.design_3stage(soft_iters=50, temp_iters=25, hard_iters=10)

            name = "pmhc_%s_%03d" % (arm, n)
            pdb_out = "%s/%s.pdb" % (W, name)
            model.save_pdb(pdb_out)

            aux = model.aux
            pae = aux.get("pae")
            pae_saved = False
            if pae is not None:
                np.save("%s/%s_pae.npy" % (W, name), np.asarray(pae))
                pae_saved = True

            met = measure(pdb_out, crop_a, crop_c)
            log = aux.get("log", {})
            rec = dict(name=name, arm=arm, seed=seed, length=L,
                       sequence=model.get_seq()[0], hotspot=hot,
                       num_pos=cfg["num_pos"], pae_saved=pae_saved,
                       i_ptm=round(float(log.get("i_ptm", 0)), 3),
                       plddt=round(float(log.get("plddt", 0)), 3), **met)
            records.append(rec)
            print("%-18s len %2d | ON-MUTATION %d (%d atoms) | peptide %2d groove %2d "
                  "| i_ptm %.2f plddt %.2f | pae %s"
                  % (name, L, rec["on_mutation"], rec["mutation_contacts"],
                     rec["on_peptide"], rec["on_groove"], rec["i_ptm"], rec["plddt"],
                     "saved" if pae_saved else "MISSING"))
            n += 1

    json.dump(records, open("%s/pmhc_probe.json" % W, "w"), indent=0)
    print("\nwrote %s/pmhc_probe.json  (%d designs, %.1f min)"
          % (W, len(records), (time.time() - t0) / 60.0))

    print("\n==================== PROBE RESULT ====================")
    print("THE QUESTION: does the tight set (5,6,7) reach residue 6 more often than")
    print("the obvious set (4-9)? Phase 1 argued yes. This measures it.\n")
    print("%-10s %5s %12s %10s %10s %9s"
          % ("arm", "n", "reached mut", "mean on-mut", "mean pept", "mut share"))
    for arm in ARMS:
        rs = [r for r in records if r["arm"] == arm]
        if not rs:
            print("%-10s     0   (no designs completed)" % arm)
            continue
        reached = sum(1 for r in rs if r["on_mutation"] > 0)
        mean_mut = np.mean([r["on_mutation"] for r in rs])
        mean_pep = np.mean([r["on_peptide"] for r in rs])
        share = np.mean([r["on_mutation"] / r["on_peptide"]
                         for r in rs if r["on_peptide"] > 0]) if any(
                             r["on_peptide"] > 0 for r in rs) else float("nan")
        print("%-10s %5d %7d/%-4d %10.2f %10.2f %8.0f%%"
              % (arm, len(rs), reached, len(rs), mean_mut, mean_pep, 100 * share))
    print("\nPAE saved for %d of %d designs. ipSAE is computable for the first time"
          % (sum(1 for r in records if r["pae_saved"]), len(records)))
    print("in this project if that number is not zero.")

    print("\n===== COPY FROM HERE =====")
    for r in records:
        print(json.dumps(r))
    print("===== COPY TO HERE =====")


main()
