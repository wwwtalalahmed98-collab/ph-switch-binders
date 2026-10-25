# Our ECD epitope against the one Adaptyv hands out

**25 Oct 2026.** Adaptyv published `design-a-protein.com`, a guided on-ramp for this
challenge. It offers one GLP-1R target and it is **3IOL** — the isolated extracellular
domain with GLP-1 bound — with the user picking hotspots on that surface and
RFD4-Proteina / Protpardelle-1c + ProteinMPNN generating a binder.

We designed against **6X18**: the full-length receptor, cryo-EM, GLP-1 bound. Different
structure, different method, different year. That makes the comparison worth running:
if the two disagree about which receptor residues the hormone actually grips, we are
designing against a surface the lab assay may not present.

## The comparison

Receptor residues within 4.5 Å of the hormone, UniProt numbering, ECD only (≤145).

| | residues |
|---|---|
| 6X18, ours | 22 |
| 3IOL, theirs | 17 |
| shared | **16** |

```
shared (16)   V30 S31 L32 T35 V36 W39 D67 E68 Y69 Y88 L89 P90 W91 R121 L123 E128
6X18 only (6) R43 P137 E138 L141 L144 Y145
3IOL only (1) T29
```

- **16 of their 17 contacts (94%) are in our set.**
- 16 of our 22 (73%) are in theirs.
- Jaccard overlap 0.70.
- **Zero residue-identity mismatches** at the 16 shared positions — the amino acid at
  each number is the same in both files.

## What it means, and what it does not

The agreement is the useful part: our landing-pad binders are aimed at the groove the
official on-ramp points at. That is independent of our own reasoning, which is what
makes it worth anything.

The six extras are not noise. P137, E138, L141, L144 and Y145 are the ECD's
C-terminal stalk, which runs down toward the transmembrane bundle and **only exists in
a full-length structure**. 3IOL is cut off before it. So our designs can engage a
stretch of receptor that anyone designing from 3IOL cannot see — relevant because that
stalk is on the path a two-domain molecule has to travel.

The zero-mismatch result is a numbering check, and this project has had three numbering
bugs. Two independently deposited structures agreeing residue-for-residue on 16
positions means our UniProt mapping is right.

**What this does not establish:** nothing here says our binders bind, and nothing here
says they activate. It says they are pointed at the correct piece of protein. The lab
decides the rest.

Reproduce: `python glp1r/crosscheck_3iol.py` (downloads nothing; expects
`data/glp1r/6X18.pdb` and `data/glp1r/3IOL.pdb`).
