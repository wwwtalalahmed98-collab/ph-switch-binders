# Running BindCraft for Challenge 1 — step by step

**Notebook:** https://colab.research.google.com/github/martinpacesa/BindCraft/blob/main/notebooks/BindCraft.ipynb

## Files to upload

From `targets/`:

| File | Purpose |
|---|---|
| `EGFR_A_domainI_ligand.pdb` | Epitope A target — domain I, 161 residues (**run this first**) |
| `EGFR_B_domainII_dimerarm.pdb` | Epitope B target — domain II, 146 residues (only if time allows) |
| `bindcraft_EGFR_A.json` | Target settings for epitope A |
| `bindcraft_EGFR_B.json` | Target settings for epitope B |

Upload the `.pdb` files to `/content/` in the Colab session (drag into the file browser,
or use the upload cell). Mount Google Drive so results survive a disconnect.

## Settings

**Epitope A (primary):**
- `starting_pdb`: `/content/EGFR_A_domainI_ligand.pdb`
- `chains`: `A`
- `target_hotspot_residues`: `20,21,22,25,28,50,51`
- `lengths`: `[55, 95]`
- `number_of_final_designs`: `12`

**Epitope B (hedge):**
- `target_hotspot_residues`: `250,253,254,255,257,258,259`

**Presets:**
- Advanced: `default_4stage_multimer.json`
- Filters: **`relaxed_filters.json`** (not `default_filters.json` — see throughput note)

## Why these numbers

- **Hotspots** are the Phase 1 epitope, in PDB numbering (UniProt − 24). Epitope A is
  F44, E45, D46, L49, Q52, Y74, D75 — 100% identical human↔mouse, overlapping the EGF
  ligand site, and verified to have identical solvent exposure in the trimmed construct
  and the full ECD (no occlusion, no trimming artefact).
- **Y275 is deliberately excluded** from the epitope B hotspots. It gains 0.30 → 0.63
  relative SASA on trimming, meaning it is half-occluded by domain IV in the intact
  tethered receptor. Targeting it would produce binders that clash with domain IV on the
  protein Adaptyv actually assays.
- **Targets are trimmed to a single domain** (161 and 146 residues). BindCraft's own
  documentation: *"Always try to trim the input target PDB to the smallest size possible!
  It will significantly speed up the binder generation and minimise the GPU memory
  requirements."* The full ECD is 609 residues and will not run on a free GPU.
- **Lengths 55–95** keep every design inside the competition's **minibinder** category
  (40–100 aa, FAQ 2), which has the best tool support and expression rates, and avoids
  straddling two scoring categories.

## Throughput reality — read this before starting

BindCraft's README states Colab is *not recommended* and that you should expect "at least
a few hundred trajectories to see some accepted binders," with 100 filtered designs as the
normal target before picking 5–20 to order. **That is not achievable on a free T4 by
Oct 4.** The adaptations:

1. **Trimmed targets** (done) — the single biggest speedup available.
2. **`relaxed_filters.json`** instead of default filters, so trajectories are not discarded
   at a bar calibrated for labs generating hundreds of designs. We re-rank ourselves in
   Phase 4 using the organizers' own metric, so the filter is not our last line of defence.
3. **Accept fewer designs.** We need 20 total, not 100. Anything from ~12–15 accepted
   designs at epitope A is enough to fill the submission after the controls.

**Strongly consider Colab Pro (~$10).** It buys longer runtimes and an L4/A100 instead of a
T4, which is likely the difference between ~5 designs and ~20 by the deadline. Given the
deadline is Oct 4, this is the highest-leverage $10 in the whole project. Modal's free tier
is the alternative, but the container setup cost eats hours we don't have.

## While it runs

Save everything to Drive. Send me the output CSV / accepted design folder as it fills up —
I can start Phase 3 (histidine placement against the conserved carboxylates E45, D46, D75,
E97) on the first accepted designs without waiting for the full set.

## If BindCraft stalls out

Fallback is RFdiffusion + ProteinMPNN + AF2 validation, which is often cheaper per design
on a weak GPU because it is a single forward pass per candidate rather than iterative
hallucination. Say the word and I'll prepare those inputs from the same target constructs.
