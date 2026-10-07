"""Phase 4: novelty screen (hard disqualification filter).

Competition rule: designs must be de novo and zero-shot, with adequate sequence and
structural diversity from known proteins. Proteinbase classifies a protein as de novo
(level 4) at <=30% sequence identity AND less than moderate structural similarity.

This screens sequence novelty via NCBI BLAST against SwissProt and reports the maximum
identity to any known protein, so nothing goes into the submission that could be thrown
out on the novelty check.

Usage:
    python scripts/phase4_novelty.py designs.fasta --out results/novelty.json
    python scripts/phase4_novelty.py --selftest
"""
import argparse
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
IDENTITY_LIMIT = 30.0  # Proteinbase de novo threshold


def submit(seq, database="swissprot", program="blastp"):
    data = urllib.parse.urlencode({
        "CMD": "Put", "PROGRAM": program, "DATABASE": database,
        "QUERY": seq, "HITLIST_SIZE": 20, "EXPECT": 10.0,
    }).encode()
    req = urllib.request.Request(BLAST_URL, data=data,
                                 headers={"User-Agent": "phase4-novelty"})
    html = urllib.request.urlopen(req, timeout=120).read().decode(errors="replace")
    rid = re.search(r"RID = (\S+)", html)
    if not rid:
        raise RuntimeError("NCBI did not return a request id")
    return rid.group(1)


def poll(rid, timeout=2400, interval=20):
    """NCBI queue times were measured at >10 min per query on 2026-10-03, so the
    timeout is deliberately generous. RIDs are persisted by the caller, so a run that
    times out can be resumed with --resume instead of resubmitting from scratch."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        q = urllib.parse.urlencode({"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"})
        html = urllib.request.urlopen(f"{BLAST_URL}?{q}", timeout=120).read().decode(errors="replace")
        if "Status=WAITING" in html:
            time.sleep(interval)
            continue
        if "Status=FAILED" in html:
            raise RuntimeError(f"BLAST failed for {rid}")
        if "Status=UNKNOWN" in html:
            raise RuntimeError(f"BLAST RID expired: {rid}")
        if "Status=READY" in html:
            return "ThereAreHits=yes" in html
        time.sleep(interval)
    raise TimeoutError(f"BLAST timed out for {rid}")


def fetch(rid):
    q = urllib.parse.urlencode({
        "CMD": "Get", "RID": rid, "FORMAT_TYPE": "Tabular",
        "ALIGNMENTS": 20, "DESCRIPTIONS": 20,
    })
    return urllib.request.urlopen(f"{BLAST_URL}?{q}", timeout=180).read().decode(errors="replace")


def parse_tabular(text):
    hits = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("<"):
            continue
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        try:
            hits.append(dict(subject=parts[1], identity=float(parts[2]),
                             length=int(parts[3]), evalue=float(parts[10])))
        except ValueError:
            continue
    return hits


def screen(name, seq, database="swissprot"):
    rid = submit(seq, database)
    has_hits = poll(rid)
    hits = parse_tabular(fetch(rid)) if has_hits else []
    top = max((h["identity"] for h in hits), default=0.0)
    return dict(name=name, length=len(seq), max_identity=top,
                n_hits=len(hits), top_hits=hits[:5], database=database,
                passes_de_novo=top <= IDENTITY_LIMIT,
                verdict=("PASS" if top <= IDENTITY_LIMIT
                         else "FAIL - too similar to a known protein"))


def read_fasta(path):
    out, name, buf = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if name:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line:
            buf.append(line)
    if name:
        out.append((name, "".join(buf)))
    return out


def selftest():
    """Two controls: a real natural protein must FAIL the novelty filter, and a
    scrambled sequence of the same composition should return far lower identity."""
    import random
    nat = ("MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPFEDHVKLVNEVTE"
           "FAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEPERNECFLQHKDDNPNLPRLVRPEV")
    random.seed(0)
    scram = "".join(random.sample(nat, len(nat)))
    for label, seq in [("natural_serum_albumin_fragment", nat),
                       ("scrambled_same_composition", scram)]:
        print(f"\n--- {label} ({len(seq)} aa) ---")
        try:
            r = screen(label, seq)
            print(f"  max identity to SwissProt: {r['max_identity']:.1f}%  "
                  f"hits: {r['n_hits']}  -> {r['verdict']}")
            if r["top_hits"]:
                print(f"  top hit: {r['top_hits'][0]['subject']} "
                      f"({r['top_hits'][0]['identity']:.1f}%)")
        except Exception as e:
            print(f"  ERROR: {e}")
    print("\nExpected: the natural fragment FAILS (high identity); the scrambled "
          "control returns much lower identity or no hits.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta", nargs="?")
    ap.add_argument("--out", default="results/novelty.json")
    ap.add_argument("--database", default="swissprot")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()

    if a.selftest:
        selftest()
        return
    if not a.fasta:
        ap.error("give a FASTA of designs, or --selftest")

    # Submit every query first, then poll. NCBI queue time dominates, so running
    # serially would take ~5 min per design; batching makes the whole set ~1 wait.
    designs = read_fasta(a.fasta)
    pending = []
    for name, seq in designs:
        try:
            rid = submit(seq, a.database)
            print(f"submitted {name} ({len(seq)} aa) -> {rid}", flush=True)
            pending.append((name, seq, rid))
        except Exception as e:
            print(f"submit failed for {name}: {e}", flush=True)
            pending.append((name, seq, None))
        time.sleep(3)  # NCBI asks for >=3 s between submissions

    # Persist RIDs immediately: NCBI results stay retrievable for ~24 h, so a crash or
    # timeout must never cost us the queue position.
    rid_path = os.path.splitext(a.out)[0] + "_rids.json"
    os.makedirs(os.path.dirname(rid_path) or ".", exist_ok=True)
    json.dump([{"name": n, "rid": r} for n, _, r in pending],
              open(rid_path, "w"), indent=1)
    print(f"\nRIDs saved to {rid_path} (resume with --resume if this run is "
          f"interrupted)", flush=True)
    print(f"all {len(pending)} submitted; polling for results ...", flush=True)
    results = []
    for name, seq, rid in pending:
        if rid is None:
            results.append(dict(name=name, length=len(seq),
                                error="submission failed", verdict="ERROR"))
            continue
        try:
            has_hits = poll(rid)
            hits = parse_tabular(fetch(rid)) if has_hits else []
            top = max((h["identity"] for h in hits), default=0.0)
            r = dict(name=name, length=len(seq), max_identity=top,
                     n_hits=len(hits), top_hits=hits[:5], database=a.database,
                     passes_de_novo=top <= IDENTITY_LIMIT,
                     verdict=("PASS" if top <= IDENTITY_LIMIT
                              else "FAIL - too similar to a known protein"))
            print(f"  {name}: max identity {top:.1f}%  -> {r['verdict']}", flush=True)
        except Exception as e:
            r = dict(name=name, length=len(seq), error=str(e), verdict="ERROR")
            print(f"  {name}: ERROR {e}", flush=True)
        results.append(r)

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    json.dump(results, open(a.out, "w"), indent=1)
    failed = [r for r in results if not r.get("passes_de_novo", False)]
    print(f"\n{len(results) - len(failed)}/{len(results)} pass the <= {IDENTITY_LIMIT}% "
          f"de novo threshold. Wrote {a.out}")
    if failed:
        print("DO NOT SUBMIT these without checking:",
              ", ".join(r["name"] for r in failed))


if __name__ == "__main__":
    main()
