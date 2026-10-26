"""Rank candidate epitopes on human TfR1 by human/mouse conservation.

Runs after conservation_map.py, which wrote per-residue identity and exposure.

The binding sites we must avoid or compare against are read off real complexes rather
than remembered from the literature:
  1SUV  transferrin-receptor complex      -> the transferrin site, which we avoid
  1DE4  HFE-receptor complex              -> overlaps the transferrin site, also avoided
  11RN  parvovirus VP1u on the apical tip -> the recommended apical epitope, for comparison

A patch is every exposed residue whose CA lies within PATCH_RADIUS of a centre residue.
It is scored by the fraction of its exposed residues identical between human and mouse.
That is the quantity Challenge 4 objective 2 depends on, and nothing else here is a proxy
for it.
"""
import io
import json
import os

import numpy as np
import requests
from Bio import PDB

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "data", "tfr1")
PATCH_RADIUS = 10.0
CONTACT = 5.0
MIN_PATCH = 8          # a patch smaller than this is not an interface
AVOID_PAD = 8.0        # keep this far from the transferrin / HFE footprint

THREE2ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
             "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
             "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
             "TYR": "Y", "VAL": "V"}


def cif(pid):
    p = os.path.join(D, pid + ".cif")
    if not os.path.exists(p):
        r = requests.get("https://files.rcsb.org/download/%s.cif" % pid, timeout=60)
        r.raise_for_status()
        io.open(p, "w", encoding="utf-8", newline="\n").write(r.text)
    return PDB.MMCIFParser(QUIET=True).get_structure(pid, p)[0]


def receptor_chain_ids(model, human):
    """Which chains are TfR1? Decided by matching the actual residues to the human
    sequence, not by assuming a chain letter."""
    out = []
    for c in model:
        hit = tot = 0
        for r in c:
            if r.id[0] != " " or r.get_resname() not in THREE2ONE:
                continue
            n = r.id[1]
            if 1 <= n <= len(human):
                tot += 1
                hit += THREE2ONE[r.get_resname()] == human[n - 1]
        if tot > 80 and hit / float(tot) > 0.9:
            out.append(c.id)
    return out


def partner_contacts(pid, human):
    """TfR1 residue numbers contacted by whatever else is in the structure."""
    m = cif(pid)
    rec = set(receptor_chain_ids(m, human))
    if not rec:
        return set(), []
    ra, pa = [], []
    for c in m:
        for r in c:
            if r.id[0] != " ":
                continue
            for a in r:
                if a.element == "H":
                    continue
                (ra if c.id in rec else pa).append((a.get_coord(), r.id[1], c.id))
    if not pa:
        return set(), sorted(rec)
    rc = np.array([x[0] for x in ra])
    pc = np.array([x[0] for x in pa])
    hit = set()
    for i in range(0, len(pc), 500):
        d = np.linalg.norm(rc[:, None, :] - pc[None, i:i + 500, :], axis=2)
        for j in np.where((d < CONTACT).any(axis=1))[0]:
            hit.add(ra[j][1])
    return hit, sorted(rec)


def main():
    human = "".join(l.strip() for l in io.open(os.path.join(D, "P02786.fasta"))
                    if not l.startswith(">"))
    cons = json.load(open(os.path.join(D, "tfr1_conservation.json")))
    res = {r["uniprot"]: r for r in cons["residues"]}
    print("ectodomain identity %.1f%%" % (100 * cons["ectodomain_identity"]))

    sites = {}
    for pid in ("1SUV", "1DE4", "11RN"):
        hits, rec = partner_contacts(pid, human)
        sites[pid] = hits
        print("%s: receptor chains %s, %d contacted residues" % (pid, rec, len(hits)))
        if hits:
            print("   %s" % ", ".join("%s%d" % (res[n]["aa"], n) for n in sorted(hits)
                                      if n in res)[:300])

    apical = sites["11RN"]
    if apical:
        got = [n for n in apical if n in res]
        c = sum(res[n]["conserved"] for n in got)
        print("\nRECOMMENDED APICAL EPITOPE (from 11RN): %d of %d residues conserved"
              " human/mouse = %.0f%%" % (c, len(got), 100.0 * c / max(len(got), 1)))
        print("   differing: %s" % ", ".join(
            "%s%d" % (res[n]["aa"], n) for n in sorted(got) if not res[n]["conserved"]))

    avoid_xyz = np.array([res[n]["xyz"] for n in (sites["1SUV"] | sites["1DE4"])
                          if n in res]) if (sites["1SUV"] | sites["1DE4"]) else None

    exp = [r for r in res.values() if r["rel_sasa"] >= 0.20]
    xyz = np.array([r["xyz"] for r in exp])
    rows = []
    for i, centre in enumerate(exp):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        mem = [exp[j] for j in np.where(d <= PATCH_RADIUS)[0]]
        if len(mem) < MIN_PATCH:
            continue
        if avoid_xyz is not None:
            if np.linalg.norm(avoid_xyz - xyz[i], axis=1).min() < AVOID_PAD:
                continue
        nc = sum(m["conserved"] for m in mem)
        rows.append(dict(centre=centre["uniprot"], centre_aa=centre["aa"],
                         n=len(mem), conserved=nc, frac=nc / float(len(mem)),
                         mean_sasa=round(float(np.mean([m["rel_sasa"] for m in mem])), 3),
                         members=sorted(m["uniprot"] for m in mem),
                         n_his=sum(m["aa"] == "H" for m in mem),
                         n_basic=sum(m["aa"] in "KR" for m in mem),
                         n_acid=sum(m["aa"] in "DE" for m in mem)))

    rows.sort(key=lambda r: (r["frac"], r["n"]), reverse=True)
    print("\n%d candidate patches clear of the transferrin/HFE footprint" % len(rows))
    print("\n%-9s %4s %5s %6s %5s %4s %4s  %s"
          % ("centre", "size", "cons", "frac", "sasa", "His", "K/R", "span"))
    for r in rows[:18]:
        print("%-9s %4d %5d %5.0f%% %5.2f %4d %4d  %d-%d"
              % ("%s%d" % (r["centre_aa"], r["centre"]), r["n"], r["conserved"],
                 100 * r["frac"], r["mean_sasa"], r["n_his"], r["n_basic"],
                 min(r["members"]), max(r["members"])))

    out = os.path.join(D, "tfr1_patches.json")
    json.dump(dict(patches=rows[:60], apical=sorted(apical),
                   transferrin_site=sorted(sites["1SUV"] | sites["1DE4"])),
              open(out, "w"), indent=0)
    print("\nwrote %s" % out)


if __name__ == "__main__":
    main()
