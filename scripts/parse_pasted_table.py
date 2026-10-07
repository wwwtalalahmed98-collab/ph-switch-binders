"""Parse the tab-separated table printed by the in-Kaggle switch analysis cell.

Accepts the block between the COPY markers (markers may be included or stripped) and
emits designs.json in the schema phase5_submission.py expects. Tolerates the mangling
that copy-paste through a chat client tends to introduce: lost tabs collapsed to runs of
spaces, blank lines, and stray wrapping.

Usage:
    python scripts/parse_pasted_table.py pasted.txt --out results/designs.json
"""
import argparse
import json
import os
import re
import sys

COLUMNS = ["name", "role", "n_his_contacts", "switch_score",
           "his_substitutions", "plddt", "sc_rmsd", "sequence"]
AA = set("ACDEFGHIKLMNPQRSTVWY")


def backbone_of(name):
    """egfr_a_0_s3_sw -> egfr_a_0.

    Strip the role suffix first: a plain rsplit on '_s' would cut at the '_sw' suffix
    and leave the sequence index attached, making every sequence look like its own
    backbone and silently defeating the diversity cap in phase5_submission.
    """
    for suffix in ("_sw", "_ctl"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return re.sub(r"_s\d+$", "", name)


def num(v, cast=float):
    v = (v or "").strip()
    if v in ("", "None", "null", "NA", "-"):
        return None
    try:
        return cast(v)
    except ValueError:
        return None


def split_row(line):
    """Prefer real tabs; fall back to runs of 2+ spaces if the paste lost them."""
    if "\t" in line:
        parts = [p.strip() for p in line.split("\t")]
    else:
        parts = [p.strip() for p in re.split(r"\s{2,}", line.strip())]
    if len(parts) < len(COLUMNS):
        # last resort: single spaces, but only if the field count then matches
        alt = [p.strip() for p in line.split()]
        if len(alt) == len(COLUMNS):
            parts = alt
    return parts


def parse(text):
    rows, skipped = [], []
    for raw in text.splitlines():
        line = raw.rstrip()
        if not line.strip():
            continue
        if "COPY FROM HERE" in line or "COPY TO HERE" in line:
            continue
        if line.lower().startswith("name\t") or line.split()[:1] == ["name"]:
            continue  # header
        parts = split_row(line)
        if len(parts) < len(COLUMNS):
            skipped.append((raw, f"{len(parts)} fields, expected {len(COLUMNS)}"))
            continue
        # sequence is the last field; anything extra belongs to substitutions
        rec = dict(zip(COLUMNS[:-1], parts[:len(COLUMNS) - 1]))
        rec["sequence"] = parts[-1].strip().upper()
        bad = set(rec["sequence"]) - AA
        if bad or not rec["sequence"]:
            skipped.append((raw, f"bad residues {''.join(sorted(bad)) or 'empty'}"))
            continue
        rows.append(dict(
            name=rec["name"],
            sequence=rec["sequence"],
            role=rec["role"],
            n_his_contacts=num(rec["n_his_contacts"], int) or 0,
            switch_score=num(rec["switch_score"]) or 0.0,
            his_substitutions=rec["his_substitutions"],
            plddt=num(rec["plddt"]),
            sc_rmsd=num(rec["sc_rmsd"]),
            epitope_conservation=1.0,
            generator="RFdiffusion",
            source_backbone=backbone_of(rec["name"]),
        ))
    return rows, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", nargs="?", help="file with the pasted table ('-' for stdin)")
    ap.add_argument("--out", default="results/designs.json")
    a = ap.parse_args()

    text = sys.stdin.read() if (not a.path or a.path == "-") else open(a.path).read()
    rows, skipped = parse(text)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(rows, open(a.out, "w"), indent=1)

    sw = [r for r in rows if r["role"] == "switch"]
    ct = [r for r in rows if r["role"] == "control"]
    print(f"parsed {len(rows)} designs -> {a.out}")
    print(f"   {len(sw)} switch, {len(ct)} control")
    if rows:
        lens = [len(r["sequence"]) for r in rows]
        print(f"   length {min(lens)}-{max(lens)} aa")
        bbs = {r["source_backbone"] for r in rows}
        print(f"   {len(bbs)} distinct backbones")
        with_q = [r for r in rows if r["plddt"] is not None]
        if with_q:
            print(f"   pLDDT {min(r['plddt'] for r in with_q):.1f}-"
                  f"{max(r['plddt'] for r in with_q):.1f}")
        hc = [r["n_his_contacts"] for r in sw]
        if hc:
            print(f"   His contacts on switch designs: "
                  f"{min(hc)}-{max(hc)}")
    fa = os.path.splitext(a.out)[0] + "_for_blast.fasta"
    with open(fa, "w") as fh:
        for r in rows:
            fh.write(f">{r['name']}\n{r['sequence']}\n")
    print(f"   wrote {fa}")
    if skipped:
        print(f"\n{len(skipped)} line(s) skipped:")
        for raw, why in skipped[:5]:
            print(f"   [{why}] {raw[:80]}")


if __name__ == "__main__":
    main()
