# Does the designed interface survive in the full extracellular region?

Adaptyv measures affinity against the **full** EGFR extracellular region. Our designs
were generated against a trimmed domain I construct (UniProt 28-188). Trimming a
domain out of a larger protein exposes surface that is buried in the intact
molecule, and a binder is perfectly happy to exploit that surface - it simply will
not be there in the assay.

Comparing rel_SASA residue by residue, isolated construct vs chain A of 6ARU:

- 161 domain I residues compared
- **10** are made artificially exposed by the trim (delta rel_SASA > 0.25):
  UniProt 105, 138, 141, 142, 158, 161, 162, 163, 166, 188

Against the contact spans the binders actually used:

| span | residues | trim-only exposure | buried in full ECD, exposed in trim |
|---|---|---|---|
| 38-152 (most designs) | 115 | 4 | 4 - (108, 136, 141, 142) |
| 38-184 (r2_w1_004) | 147 | 9 | 6 - (108, 136, 141, 142, 162, 166) |
| 28-62 (r2_w0_005) | 35 | 0 | 0 |
| 69-114 (r2_w0_004) | 46 | 1 | 1 - (108) |

So the risk is real but contained: the artefact residues are a small minority of each
span, and the span is itself much larger than the ~19-27 residue contact set. What
this does **not** yet establish is whether any individual design leans on those
specific residues, which requires its own contact list rather than the span. That
check runs over the per-design `epitope_residues` before anything is ranked, and any
design whose interface depends materially on 105/136/138/141/142/158/161/162/163/166
or 108 should be demoted regardless of its i_ptm.

Phase 2 had already shown the *epitope A hotspots* keep identical rel_SASA on
trimming. That was a check of the epitope we intended to use, and the binders did not
use it, so it does not transfer. This check is on the surface they actually chose.
