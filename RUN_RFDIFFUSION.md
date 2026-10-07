# RFdiffusion fallback — step by step

Use this if BindCraft is too slow on the available GPU. It is often cheaper per design on a
weak GPU: diffusion is a fixed number of steps per backbone, then one ProteinMPNN pass and
one AlphaFold check, rather than BindCraft's iterative hallucination loop.

**Notebook:** https://colab.research.google.com/github/sokrypton/ColabDesign/blob/main/rf/examples/diffusion.ipynb

Same target constructs as BindCraft — upload `targets/EGFR_A_domainI_ligand.pdb`
(and `EGFR_B_domainII_dimerarm.pdb` if you get to epitope B).

## Cell 3 — run RFdiffusion

**Epitope A (primary):**

| Field | Value |
|---|---|
| `name` | `EGFRa` |
| `contigs` | `A1-161:55-95` |
| `pdb` | `EGFR_A_domainI_ligand.pdb` (or leave blank for the upload prompt) |
| `hotspot` | `A20,A21,A22,A25,A28,A50,A51` |
| `iterations` | `50` |
| `num_designs` | `32` |
| `use_beta_model` | run **two batches** — see below |

**Epitope B (hedge):**

| Field | Value |
|---|---|
| `name` | `EGFRb` |
| `contigs` | `A165-310:55-95` |
| `hotspot` | `A250,A253,A254,A255,A257,A258,A259` |
| `num_designs` | `16` |

Syntax note: `A1-161:55-95` means "keep target chain A residues 1–161, diffuse a binder of
length 55–95 sampled randomly." The `55-95` range keeps every design in the competition's
**minibinder** category (40–100 aa).

### Run two batches at epitope A, one with `use_beta_model=True`

The notebook's own hint: switch to the beta parameters "if you are seeing lots of helices,
for a better SSE balance." This is worth doing deliberately, not just as a fix.

In Anthropic's campaign only **15 of 354 binders** had design models containing ≥20%
β-strand — nearly everything that works is all-α, because that is what the methods produce
by default. Track 3 selection explicitly weighs **design novelty**, so a minority of
β-containing backbones at a conserved epitope is a cheap way to be distinguishable from the
crowd. Suggested split: ~2/3 of trajectories with `use_beta_model=False`, ~1/3 with `True`.

A second reason: β-sheet interfaces present backbone and side-chain geometry that is easier
to position precisely against a target carboxylate, which is what the Phase 3 histidine
placement needs.

## Cell 5 — ProteinMPNN + AlphaFold validation

| Field | Value | Why |
|---|---|---|
| `num_seqs` | `8` | 8 sequences per backbone |
| `mpnn_sampling_temp` | `0.1` | low diversity, higher quality |
| `rm_aa` | `C` | no cysteines — avoids disulfide scrambling and expression problems |
| `use_solubleMPNN` | **`True`** | Anthropic used SolubleMPNN for 1,133 of 1,315 tested designs |
| `initial_guess` | **`True`** | the notebook recommends this for binder design |
| `num_recycles` | `3` | recommended alongside initial_guess |
| `use_multimer` | `True` | AF2-Multimer v3 params |

**Filter:** keep designs with **pae_interaction < 10** (the standard cutoff from the AF2
initial-guess binder paper) and pLDDT > 80. If almost nothing passes, keep the best by
pae_interaction and record that the cutoff was relaxed — do not silently move the bar.

Do **not** let ProteinMPNN place the switch histidines. Sequences come out of this stage
pH-naive; histidines go in during Phase 3, positioned against the specific conserved
carboxylates (E45, D46, D75, E97 at epitope A), where their geometry can be checked.

## If you get a real GPU (Modal, Colab Pro A100, cluster)

Official RFdiffusion CLI, same parameters:

```bash
python scripts/run_inference.py inference.output_prefix=outputs/EGFRa inference.input_pdb=EGFR_A_domainI_ligand.pdb 'contigmap.contigs=[A1-161/0 55-95]' 'ppi.hotspot_res=[A20,A21,A22,A25,A28,A50,A51]' inference.num_designs=100 denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0
```

For the β-balanced batch add:

```bash
inference.ckpt_override_path=models/Complex_beta_ckpt.pt
```

`noise_scale_ca=0 noise_scale_frame=0` is the standard setting for binder design — noise-free
sampling measurably improves success rates for PPI targets.

## What to send me

The accepted designs plus their scores (pae_interaction, pLDDT, ipTM). Sequences alone are
enough to start Phase 3; the PDBs let me check histidine geometry against the target
carboxylates properly, so send those if they are easy to export.

## Which pipeline to prefer

Run BindCraft first — its hit rates are better and it needs no post-filtering judgement.
Switch here if you are getting fewer than roughly one accepted design per hour, or if
BindCraft runs out of memory. Running both and pooling is also fine, and is actively
useful: the organizers' own protocol required designs to be drawn from at least three
generators with no single method contributing more than half, precisely to avoid
method-specific failure modes.
