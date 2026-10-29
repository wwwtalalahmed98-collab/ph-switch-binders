"""Re-screen the MPNN pass-2 candidates with a proline/glycine gate.

The two-sided screen in cellTFR1j_mpnn2.py checks charge, hydrophobicity,
Kyte-Doolittle and identical-residue runs. It does NOT look at proline, so a
sequence that is 30% proline - a polyproline rod, not a globular fold - passes
as "good". Natural globular proteins are about 5% proline and 7% glycine.

This script also recomputes the composition numbers the Kaggle run printed. A
mismatch means the candidate table was transcribed wrong, which matters because
that table was copied out of a log by hand.
"""
import io, sys

KD = {"A":1.8,"R":-4.5,"N":-3.5,"D":-3.5,"C":2.5,"Q":-3.5,"E":-3.5,"G":-0.4,
      "H":-3.2,"I":4.5,"L":3.8,"K":-3.9,"M":1.9,"F":2.8,"P":-1.6,"S":-0.8,
      "T":-0.7,"W":-0.9,"Y":-1.3,"V":4.2}
CHARGED = set("DEKR")
HYDRO = set("AVILMFWYC")

# proline and glycine are both fold-breakers in bulk: proline cannot donate a
# backbone hydrogen bond and kinks helices, glycine costs the chain its rigidity.
P_MAX, G_MAX, PG_MAX = 0.12, 0.12, 0.20

rows = []
with io.open("tfr1/results/mpnn2_candidates.tsv", encoding="utf-8") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))

print("%-22s %4s %5s %5s %5s %6s  %-8s %s"
      % ("name", "on", "P%", "G%", "P+G%", "kd", "was", "now"))
survivors = []
for r in rows:
    s = r["seq"]; L = len(s)
    p, g = s.count("P") / L, s.count("G") / L
    kd = sum(KD[c] for c in s) / L
    ok = p <= P_MAX and g <= G_MAX and (p + g) <= PG_MAX
    now = r["tier"] if ok else "fold-breaker"
    if ok and r["tier"] == "good":
        survivors.append(r)
    print("%-22s %4s %4.0f%% %4.0f%% %4.0f%% %+6.2f  %-8s %s"
          % (r["name"], r["on_site_parent"], 100*p, 100*g, 100*(p+g), kd,
             r["tier"], now))

print("\n%d of %d candidates survive the proline/glycine gate."
      % (sum(1 for r in rows if r["seq"].count("P")/len(r["seq"]) <= P_MAX
             and r["seq"].count("G")/len(r["seq"]) <= G_MAX
             and (r["seq"].count("P")+r["seq"].count("G"))/len(r["seq"]) <= PG_MAX),
         len(rows)))
print("%d of those were also tier 'good': %s"
      % (len(survivors), ", ".join(r["name"] for r in survivors)))
