"""Write the sealed GLP-1R target as a CIF with ONE numbering system.

Biopython's MMCIFIO writes label_asym_id='A' and label_seq_id=1..283 while keeping
auth_asym_id='R' and auth_seq_id=UniProt. RFD4-Proteina's selection grammar resolves
res_id against label_seq_id, so `sel('R/*/[30]')` against that file would have matched
nothing on the chain and, had the chain matched, would have selected UniProt residue 58.

Three residue-numbering bugs have already cost this project time, and every one of them
came from two numbering systems being silently different. Rather than work out which
column AtomWorks reads and trust the answer, this makes both columns the same: label and
auth agree on chain R and on UniProt residue numbers. Then any reader gets UniProt, and
a mistake shows up as a mismatch instead of as the wrong residues.

Verification at the bottom checks named hotspot residues against the file itself, not
against the map that produced it.
"""
import io
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "data", "glp1r", "GLP1R_sealed.cif")
DST = os.path.join(ROOT, "data", "glp1r", "GLP1R_sealed_rfd4.cif")
MAP = os.path.join(ROOT, "data", "glp1r", "glp1r_sealed_map.json")

lines = io.open(SRC, encoding="utf-8").read().split("\n")

cols, start = [], None
for i, ln in enumerate(lines):
    if ln.startswith("_atom_site."):
        cols.append(ln.strip().split(".", 1)[1])
        start = i + 1
ix = {c: k for k, c in enumerate(cols)}
for need in ("label_asym_id", "label_seq_id", "auth_seq_id", "auth_asym_id"):
    if need not in ix:
        raise SystemExit("source CIF lacks %s" % need)

out, changed = [], 0
for ln in lines:
    if not ln.startswith(("ATOM", "HETATM")):
        out.append(ln)
        continue
    f = ln.split()
    if len(f) != len(cols):
        raise SystemExit("row has %d fields, header declares %d" % (len(f), len(cols)))
    f[ix["label_asym_id"]] = f[ix["auth_asym_id"]]
    f[ix["label_seq_id"]] = f[ix["auth_seq_id"]]
    out.append(" ".join(f))
    changed += 1

io.open(DST, "w", encoding="utf-8", newline="\n").write("\n".join(out))
print("rewrote %d atom rows -> %s" % (changed, DST))

# --- verify against the written file, not against the map that produced it ---
seen = {}
chains = set()
for ln in io.open(DST, encoding="utf-8"):
    if not ln.startswith("ATOM"):
        continue
    f = ln.split()
    chains.add(f[ix["label_asym_id"]])
    if f[ix["label_seq_id"]] != f[ix["auth_seq_id"]]:
        raise SystemExit("label_seq_id and auth_seq_id still disagree on a row")
    seen[int(f[ix["label_seq_id"]])] = f[ix["label_comp_id"]]

print("label_asym_id values in file: %s" % sorted(chains))
print("residues: %d, numbering %d..%d" % (len(seen), min(seen), max(seen)))

T31 = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G",
           HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S",
           THR="T", TRP="W", TYR="Y", VAL="V")

m = json.load(open(MAP))
groove, channel = m["ecd_contacts"], m["tm_contacts"]
expect = {d["uniprot"]: d["aa"] for d in m["order"]}

bad = [n for n in groove + channel if T31.get(seen.get(n, "")) != expect.get(n)]
if bad:
    raise SystemExit("hotspot residues disagree with the map at %s" % bad)
print("all %d hotspots present and matching: %d groove + %d channel"
      % (len(groove) + len(channel), len(groove), len(channel)))

spot = [39, 121, 145, 190, 298, 391]
print("spot check: %s"
      % ", ".join("%s%d" % (T31.get(seen[n], "?"), n) for n in spot if n in seen))

sel_groove = "sel('R/*/[%s]')" % ", ".join(str(x) for x in groove)
sel_chan = "sel('R/*/[%s]')" % ", ".join(str(x) for x in channel)
sel_all = "sel('R/*/[%s]')" % ", ".join(str(x) for x in sorted(groove + channel))
io.open(os.path.join(ROOT, "data", "glp1r", "rfd4_selections.txt"), "w",
        encoding="utf-8", newline="\n").write(
    "groove (ECD, 22)\n%s\n\nchannel (TM, 30)\n%s\n\ncombined (52)\n%s\n"
    % (sel_groove, sel_chan, sel_all))
print("\nwrote selection strings to data/glp1r/rfd4_selections.txt")
