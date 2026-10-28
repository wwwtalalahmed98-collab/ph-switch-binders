"""Kaggle - TfR1 Challenge 4. Recover the epitope-contact measurement the design run lost.

RUNS ON CPU. No GPU quota is spent here. Attach the tfr1_generate notebook's output as an
input dataset and run this alone.

WHY THIS EXISTS. The generation run produced 47 designs and measured NOTHING about where
they bind: every one of the 47 carries

    measure_error: no chain carries the epitope numbering; {'A': 108, 'B': 68}

The check inside the worker looked for a chain whose residue NUMBERS contain the epitope
numbers (559, 560, ...). ColabDesign's save_pdb renumbers whatever it writes from 1, so the
cropped target came back as residues 1..108 and the test could never pass. The guard was
right to refuse rather than guess which chain was the target - that is the failure mode it
was written for - but it refused on every design, so nothing was measured.

The information is not lost. The crop kept a known, ordered list of UniProt residue numbers
(tfr1_targets.json, "kept"), and save_pdb writes them in that same order, so

    UniProt number of chain-A residue i  ==  kept[i - 1]

That is an ASSUMPTION about how save_pdb orders residues, so it is not taken on trust: the
amino-acid sequence of chain A is compared against the amino acids those kept numbers have
in the real human sequence, and any design whose chain A does not match is reported and
skipped rather than measured against a mapping that does not hold.

WHAT IT PRINTS. Per design: contacts on the chosen epitope and contacts anywhere else, the
same two numbers the generation run was supposed to produce. A high off-site count means the
design found the surface the crop exposed rather than the epitope, and is rejected on that.
"""
import glob
import json
import os
import urllib.request

import numpy as np
from Bio import PDB

W = "/kaggle/working"
CONTACT = 4.5          # A, any heavy atom to any heavy atom
T31 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def find(name):
    hits = glob.glob("/kaggle/input/**/%s" % name, recursive=True) + glob.glob("%s/%s" % (W, name))
    if not hits:
        raise SystemExit("cannot find %s - attach the tfr1_generate notebook output as an "
                         "input dataset" % name)
    return hits[0]


DESIGNS = json.load(open(find("tfr1_designs.json")))
TARGETS = json.load(open(find("tfr1_targets.json")))
print("%d designs, sites %s" % (len(DESIGNS), sorted(TARGETS)))

FA = "%s/P02786.fasta" % W
if not os.path.exists(FA):
    urllib.request.urlretrieve("https://rest.uniprot.org/uniprotkb/P02786.fasta", FA)
HUMAN = "".join(l.strip() for l in open(FA) if not l.startswith(">"))
print("human TfR1 %d aa" % len(HUMAN))

for tag, M in sorted(TARGETS.items()):
    print("%s: crop kept %d residues, epitope %d" % (tag, len(M["kept"]), len(M["epitope"])))


def residues(chain):
    """Ordered, de-duplicated standard residues, as save_pdb wrote them."""
    out, seen = [], set()
    for r in chain:
        if r.id[0] == " " and r.get_resname() in T31 and r.id[1] not in seen:
            seen.add(r.id[1])
            out.append(r)
    return out


def measure(pdb_path, tag):
    M = TARGETS[tag]
    kept, epitope = M["kept"], set(M["epitope"])
    model = PDB.PDBParser(QUIET=True).get_structure("c", pdb_path)[0]
    chains = {c.id: residues(c) for c in model}

    # the target is the chain whose length matches the crop AND whose sequence matches the
    # kept residues in the real human sequence; the binder is the other one
    want = "".join(HUMAN[n - 1] for n in kept)
    tid = None
    for cid, rs in chains.items():
        if len(rs) == len(kept) and "".join(T31[r.get_resname()] for r in rs) == want:
            tid = cid
            break
    if tid is None:
        got = {c: (len(rs), "".join(T31[r.get_resname()] for r in rs)[:20]) for c, rs in chains.items()}
        raise ValueError("no chain matches the crop sequence; expected %d residues starting "
                         "%s, saw %s" % (len(kept), want[:20], got))
    bid = [c for c in chains if c != tid]
    if len(bid) != 1:
        raise ValueError("expected one binder chain, got %s" % bid)
    bres = chains[bid[0]]

    ba = np.array([a.get_coord() for r in bres for a in r if a.element != "H"])
    on, off, touched = 0, 0, []
    for i, r in enumerate(chains[tid], start=1):
        uni = kept[i - 1]
        ra = np.array([a.get_coord() for a in r if a.element != "H"])
        if np.min(np.linalg.norm(ra[:, None, :] - ba[None, :, :], axis=-1)) <= CONTACT:
            if uni in epitope:
                on += 1
                touched.append(uni)
            else:
                off += 1
    rec = dict(on_site=on, off_site=off, n_epitope=len(epitope),
               on_fraction=round(on / float(max(on + off, 1)), 3), touched=sorted(touched),
               binder_chain=bid[0], target_chain=tid, binder_len=len(bres))
    cas = [r["CA"].get_coord() for r in bres if "CA" in r]
    if len(cas) > 4:
        d = [float(np.linalg.norm(cas[i] - cas[i + 4])) for i in range(len(cas) - 4)]
        rec["helicity"] = round(sum(1 for x in d if 5.5 <= x <= 7.0) / len(d), 3)
    return rec


out, failed = [], []
for d in DESIGNS:
    hits = (glob.glob("/kaggle/input/**/tfr1_%s/%s.pdb" % (d["site"], d["name"]), recursive=True)
            or glob.glob("/kaggle/input/**/%s.pdb" % d["name"], recursive=True))
    if not hits:
        failed.append((d["name"], "structure not found"))
        continue
    rec = dict(d)
    rec.pop("measure_error", None)
    try:
        rec.update(measure(hits[0], d["site"]))
    except Exception as e:
        failed.append((d["name"], "%s %s" % (type(e).__name__, str(e)[:160])))
        continue
    out.append(rec)

print("\nmeasured %d of %d designs" % (len(out), len(DESIGNS)))
for n, why in failed:
    print("  FAILED %s: %s" % (n, why))

json.dump(out, open("%s/tfr1_designs_measured.json" % W, "w"), indent=0)

out.sort(key=lambda r: (r["on_site"], r["true_iptm"]), reverse=True)
print("\n%-18s %4s %6s %6s %6s %5s %5s %6s %4s %4s" % ("design", "len", "i_ptm", "pLDDT",
      "kd", "on", "off", "on_fr", "hlx", "run"))
for r in out:
    print("%-18s %4d %6.3f %6.3f %6.2f %5d %5d %6.2f %4.2f %4d"
          % (r["name"], r["length"], r["true_iptm"], r["plddt"], r["hydropathy"],
             r["on_site"], r["off_site"], r["on_fraction"], r.get("helicity", 0.0),
             r["longest_run"]))

for cut in (2, 4, 6):
    print("%d of %d designs make at least %d epitope contacts"
          % (sum(1 for r in out if r["on_site"] >= cut), len(out), cut))

ok = [r for r in out if r["on_site"] >= 1]
if len(ok) > 4:
    r_ = float(np.corrcoef([x["true_iptm"] for x in ok], [x["on_site"] for x in ok])[0, 1])
    print("\nconfidence vs epitope contacts over %d designs: r = %+.3f" % (len(ok), r_))
    print("  on GLP-1R this was NEGATIVE (-0.96, -0.67) on a buried channel. If it is not"
          " negative here, that finding was about buried channels and does not transfer.")

print("\n===== COPY FROM HERE =====")
for r in out:
    print(json.dumps(r))
print("===== COPY TO HERE =====")
