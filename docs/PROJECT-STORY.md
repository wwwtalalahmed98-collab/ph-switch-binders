# Project story

A running, plain-English record of what this project is doing and why. Newest entry at the top.

| # | Date | Checkpoint |
|---|---|---|
| 2 | 2026-10-17 | The set is built and verified; one click left |
| 1 | 2026-10-17 | Measuring the mouse instead of guessing it |

---

## 2. 2026-10-17 — The set is built and verified; one click left

### Where we are

The twenty designs are built, checked, and ready, but I could not upload them — the final step needs your approval.

### What this project is for

We are entering Anthropic's protein design competition, Challenge 2, Track 3. We design a small protein that grabs a human inflammation protein (TNF-alpha) at normal body acidity, lets go when things turn slightly acidic, and does the same job on the mouse version. Up to 20 designs. The competition ranks three things in this order: does the acidity switch work, does it bind mouse as well as human, and how tightly does it bind.

### Why this step

Last checkpoint we had found that our way of ranking designs was wrong and had rebuilt it. One job was still running. This step reads that job, picks the final twenty, and checks them.

### What we actually did

**Read the last job.** It finished cleanly in 41 minutes. Of 29 parent designs offered acidity switches, **14 now yield a usable one** — up from 12.

**Caught a silent data trap.** Kaggle, the service we run on, shows you a *preview* of a result file, not the file. It cuts off at 102,400 bytes of our 113,975 and then tidily closes the brackets, so the truncated version reads as perfectly valid and quietly drops 3 of 29 records and corrupts a 4th. I checked those four against the job's own printed table: all four have nothing we would have used, so this set is unaffected. But reading it without checking would have looked fine and been wrong.

**Picked the ten, by a stated rule.** In the competition's own order of priorities: prefer two histidines over one (one leaves half the molecules still fully binding in acid, two leaves a quarter); prefer histidines that face a spot present in *both* species, since one facing a human-only spot cannot push off the mouse target; then prefer the better of the worse species. Six of the ten carry two histidines, up from five.

**Verified every pair by machine, not by eye.** Each switch must differ from its histidine-free twin at the histidine positions **and nowhere else** — otherwise a pH response in the lab cannot be pinned on the histidines. All 35 candidate placements passed.

**Wrote the methodology reviewers will read.** This is not a formality: a Claude-based workflow reads these write-ups when choosing which designs physically get made. The honest centre of ours is the negative result — our conservation shortcut did not work, and here is the one lineage that proves it.

### The bad news

**The job's main purpose failed.** I built it to rescue `tnf3_w1_019`, a strong design that was losing a third of its binding to a cosmetic cleanup step. I had assumed the cost came from substitutions forced by a repetitive-sequence rule, and built a cheaper, reversible version of that rule. The assumption was wrong: the damage came from replacing one lone cysteine, which is not optional. The design is out, and so is a second one for the same reason. The fix itself is correct and helped seven other lineages — but it did not do the thing it was built to do.

**Two of the ten controls still carry a native histidine of their own.** For those two pairs, a lab pH response cannot be cleanly credited to the histidine we added.

**We traded measured binding for an unmeasured argument.** Choosing two histidines over one costs real predicted binding in several designs. The two-histidine reasoning is textbook chemistry, not something we measured. If it is wrong, those were losses.

**And still: nothing has been physically tested.** Every number is a prediction. The acidity switch in particular has never been computed by anyone here — no calculation in this project models acidity at all. It is geometry plus textbook chemistry.

### How it works

Keep the earlier picture: the binder is a hand gripping a ball, and the switch is a small spring-loaded spike on one finger — flat at normal acidity, popped out in acid, pushing the hand off.

This step was choosing which finger, and how many spikes. Two spikes push off more reliably than one, but each spike also makes the grip slightly worse even when flat, so we cannot just add spikes everywhere. And a spike only works if there is something to push *against* — we now insist that the thing it pushes against exists on the mouse ball too, not just the human one. The twin design is the same hand with no spikes, so the lab can tell whether letting go was the spikes' doing.

### The diagram

```mermaid
flowchart TD
    A[Job finished<br/>29 parents offered switches] -->|35 placements hold<br/>BOTH species| B[14 parents usable]
    B -->|prefer 2 histidines,<br/>then shared push-off spot,<br/>then worse species| C[Ten chosen]
    C -->|each switch vs its twin| D{Differs ONLY at<br/>the histidines?}
    D -->|all 35 pass| E[20 designs<br/>10 pairs, written to file]
    D -->|any fail| X[Dropped loudly]
    E -->|portal checks pass| F[Upload BLOCKED<br/>needs your approval]
    style F fill:#b45309,stroke:#78350f,color:#ffffff
    style E fill:#15803d,stroke:#14532d,color:#ffffff
```

"20 designs" (green) exists on disk and passed every check. "Upload" (amber) is where we are stopped — the file is ready, the permission is not.

### What you have now

**Things that work and have produced results:**

- [tnf/results/SWITCH5_RESULTS.md](../tnf/results/SWITCH5_RESULTS.md) — the full result of the last job, the ten chosen designs with their numbers, and a table of what this replaces in the banked submission. Open it to see which five lineages got dropped and why.
- [submission/METHODS_TNF_V3.md](../submission/METHODS_TNF_V3.md) — the write-up reviewers read. Open it to see the negative finding stated plainly.
- [tnf/build_submission_v3.py](../tnf/build_submission_v3.py) — the script that picks the ten and verifies each pair. Run it again and it reproduces the same twenty designs from the same inputs.

**Built and ready, but not sent:**

- `submission/challenge2_proposed_v3.csv` — the 20 designs. Kept out of the public repository until the competition closes, so it will not appear on GitHub.
- `submission/challenge2_metrics_v3.csv` — every number behind those 20.

**Not done:** the upload itself.

### Decisions made, and why

| We chose | Instead of | Because |
|---|---|---|
| Two histidines where both survive both species | The single histidine that binds best | One histidine leaves 50% of molecules still binding in acid; two leaves 25%. Objective 1 asks for *no detectable* binding |
| Requiring the push-off partner to exist in both species | Taking the best-binding placement | A histidine facing a human-only spot gives a human-only switch, and mouse is ranked above affinity |
| Dropping `tnf3_w1_019` | Keeping it with a cheaper repair | The cost was a lone cysteine, which cannot be left in. The cheaper repair did not help |
| Verifying pairs by script | Checking the sequences by eye | I nearly mis-read a truncated file this same session; eyes are the wrong tool for this |
| Reporting the conservation shortcut as a failure in the methods | Quietly not mentioning it | Method novelty is a selection criterion, and a clean negative result is the most useful thing we found |
| Stopping at the upload | Finding another way to submit | The block is a deliberate guard on real-world actions. Routing around it is not mine to do |

### Jargon decoder

- **TNF-alpha** — a human protein that drives inflammation; what our designs grab.
- **pH** — how acidic something is. 7.4 is normal blood; 6.0 is slightly acidic, as in inflamed tissue or inside a cell's recycling compartment.
- **Histidine** — one of the twenty amino acid building blocks. It gains a positive charge right around pH 6.0, which is what makes an acidity switch possible.
- **Acid-OFF switch** — holds on at pH 7.4, lets go at pH 6.0.
- **Control / twin** — a copy of a design with the histidines left out, so a lab result can be attributed.
- **i_ptm** — the prediction model's confidence that two proteins stick together, roughly 0 to 1. **A confidence score, not a measured strength.**
- **On-site contacts** — how many of the binder's touches land on the exact patch we aimed for.
- **Lineage / parent** — one original design and all the switch versions built from it.
- **Free cysteine** — a lone sulphur-bearing residue. A manufacturing liability, so it gets replaced.
- **Novelty score** — a check the competition portal runs automatically; designs need 3 of 4 or higher.

### What is next, and what could go wrong

**Next:** you approve the upload, or do it yourself. The portal needs the designs file, the methodology, and a tick on the terms.

**Two constraints I did not know until I opened the portal:**

- **You can only submit once every 24 hours.** With the deadline on Sunday 18 October at 23:59 AoE, that means about two attempts left in total.
- **Designs need a novelty score of 3 out of 4 or higher**, checked automatically after upload. I have not been able to run that check, because it happens on upload. If any design fails it, we would need the second attempt to fix it — which is exactly why spending the first attempt soon rather than late matters.

**What could still go wrong:** the novelty check could reject some designs, and I do not know its threshold behaviour on twin pairs that differ by one or two residues. Beyond that, everything in the previous checkpoint's risk list still stands, led by the fact that the acidity switch is reasoned rather than computed.

---

## 1. 2026-10-17 — Measuring the mouse instead of guessing it

### Where we are

We found out that the way we were picking designs was wrong, fixed it, and we are one run away from replacing what we have already entered into the competition.

### What this project is for

We are entering Anthropic's protein design competition, Challenge 2, Track 3. The task is to design a small protein that sticks to a human inflammation protein called TNF-alpha at normal body acidity (pH 7.4), lets go of it when things get slightly more acidic (pH 6.0), and does the same job on the mouse version of the target so it can be tested in animals. We may submit up to 20 designs. They get made in a lab and physically tested, so the numbers we produce are predictions that will be checked against reality.

The competition ranks three things in this order: does the pH switch work, does it bind mouse as well as human, and how tightly does it bind. Mouse matters more than strength.

### Why this step

We had already submitted 20 designs. They were ranked partly on a shortcut: instead of checking whether each design actually binds mouse TNF-alpha, we counted how many of the human contact points are spelled the same way in mouse (*sequence conservation*). The reasoning was that if mouse and human look similar where the binder touches, the binder should work on both.

If that shortcut were sound, skipping the real check would save a lot of computer time. If it were not sound, we would have entered designs that fail the second-ranked objective, and we would not know.

### What we actually did

Three things, in order.

**We measured mouse binding directly.** We built a mouse version of the target from a published crystal structure (PDB entry `2TNF`) and predicted all 68 of our switch designs against it, using exactly the same model and settings as the human predictions so the two numbers are comparable.

**The shortcut turned out to carry almost no information.** Across the 11 designs that have both a shortcut score and a measured number, the two barely track each other (*Pearson correlation* r = +0.272, where 1.0 would be perfect agreement and 0 would be none). The design the shortcut ranked **first** had the **worst** measured mouse binding of the whole set.

Six of the ten switches we submitted fail the mouse objective. That includes the strongest binder in the submission, which turned out to make **zero** contacts with the mouse target.

**We found what was actually causing it.** We then tested the parent designs — the versions *before* we add the histidines that create the pH switch. 20 of 31 parents bind mouse fine. So the designs were never the problem. Comparing the 17 lineages where we have both numbers:

> The same histidine substitutions keep a median of **90% of human binding** and only **49% of mouse binding**.

That lopsidedness is not chemistry. It is our own selection rule showing through: we scored candidate histidine positions on the human structure and accepted or rejected them on human results alone. A position that human tolerates and mouse cannot was passing our filter *by construction*.

**We changed the rule and re-ran.** Now a histidine position is kept only if it holds up against both species, and each design is offered several candidate positions instead of one.

The clearest result is a single design, `tnf2_w0_007`. Same parent, same machinery, only the chosen position differs:

| histidines placed | human | mouse | mouse contacts |
|---|---|---|---|
| L30H + D27H (what the old rule picked) | 0.487 | **0.150** | **0** |
| E26H (what the new rule picks) | 0.566 | **0.678** | **15** |

A design that binds mouse well was being handed to the lab in a state where it barely touches mouse at all, because of where we put one amino acid. Three positions away is fine.

### How it works

Think of the binder as a hand gripping a ball, and the pH switch as putting a small spring-loaded spike on one finger. The spike is flat at normal acidity and pops out when things get acidic, pushing the hand off the ball. That is the switch: a histidine amino acid sitting next to a positively charged spot on the target, harmless at pH 7.4 and repulsive at pH 6.0.

The mistake was this. We tested where to put the spike by gripping only the *human* ball. Some spots on the finger are fine against the human ball but land right where the hand grips the *mouse* ball. We never checked, so we kept choosing those spots. Now we try the spike in several places and grip both balls each time, keeping only placements that work on both.

### The diagram

```mermaid
flowchart TD
    A[Designs we generated<br/>89 candidate binders] -->|31 strong enough on human| B[Parents tested on mouse]
    B -->|20 of 31 bind mouse| C[Mouse-binding lineages]
    B -.->|11 do not bind mouse| X[Set aside]
    C -->|up to 4 histidine<br/>positions tried each| D[Switch placement]
    D -->|tested against BOTH<br/>human and mouse| E{Holds both?}
    E -->|30 placements yes| F[Usable switches<br/>12 designs]
    E -->|the rest no| X
    F -->|best 10 pairs| G[Proposed submission<br/>awaiting your go-ahead]
    style D fill:#2563eb,stroke:#1e40af,color:#ffffff
    style G fill:#15803d,stroke:#14532d,color:#ffffff
```

"Switch placement" (blue) is the step we just rebuilt. "Proposed submission" (green) is where we are now — built, checked, not sent.

### What you have now

**Things that work and have produced results:**

- [tnf/results/MOUSE_RESULTS.md](../tnf/results/MOUSE_RESULTS.md) — the measured mouse numbers for all 68 designs, and the evidence that the old shortcut does not work. Open it to see the table showing six of our ten submitted switches failing.
- [tnf/results/PARENT_MOUSE_RESULTS.md](../tnf/results/PARENT_MOUSE_RESULTS.md) — the diagnosis. Shows the 90%-versus-49% finding that explains everything else.
- [tnf/results/SWITCH4_RESULTS.md](../tnf/results/SWITCH4_RESULTS.md) — results of the rebuilt placement step, including the `tnf2_w0_007` rescue.
- [tnf/results/mouse_scores.csv](../tnf/results/mouse_scores.csv) and [tnf/results/mouse_parents.csv](../tnf/results/mouse_parents.csv) — the raw numbers behind those documents, as spreadsheets.

**Built and ready, but not sent:**

- `submission/challenge2_proposed_v2.csv` — the 20 designs we would submit. Deliberately kept out of the public repository until the competition closes, so it will not appear on GitHub.
- `submission/challenge2_metrics_v2.csv` — the numbers behind each of those 20.

**Running right now:** one final job (`tnf_switch5_repair_fix`) that fixes two loose ends before we submit. See "What is next".

### Decisions made, and why

| We chose | Instead of | Because |
|---|---|---|
| Measuring mouse binding directly | Estimating it from sequence similarity | The estimate turned out to be nearly uninformative (r = +0.272), and six submitted designs fail because of it |
| Keeping a histidine only if it holds both species | Keeping it if it holds human | A position that human tolerates and mouse does not was passing the filter by construction |
| Trying several histidine positions per design | Taking the single best-scoring one | The best-scoring position is sometimes the one mouse cannot tolerate; `tnf2_w0_007` went from 0.150 to 0.678 purely by moving it |
| Preferring two histidines where both survive | Taking the single histidine the automatic rule picked | One histidine leaves half the molecules still fully active at pH 6.0; two leaves a quarter. The brief says *no detectable binding* at pH 6.0 |
| Fixing free cysteines always, breaking repetitive stretches only if it is cheap | Doing both unconditionally | One design lost a third of its binding to a cosmetic repair before the switch step even started |
| Showing you the proposed set before sending | Submitting automatically | It replaces a submission that is already banked |

### Jargon decoder

- **TNF-alpha** — a human protein that drives inflammation; the thing our designs are meant to grab.
- **Binder** — the small protein we design to stick to it.
- **pH** — how acidic something is. 7.4 is normal blood; 6.0 is slightly acidic, like inflamed or tumour tissue.
- **Histidine** — one of the twenty amino acid building blocks. Unusual in that it gains a positive charge right around pH 6.0, which is why it can act as a pH-triggered switch.
- **Acid-OFF switch** — a binder that holds on at pH 7.4 and lets go at pH 6.0.
- **i_ptm** — the confidence score the structure-prediction model gives that two proteins actually stick together. Roughly 0 to 1. **It is a confidence score, not a measured strength.**
- **AlphaFold2** — the model that predicts protein shapes and whether two proteins bind.
- **Parent / switch / control** — the original binder; the version with histidines added; and a twin of the switch with the histidines left out. The control exists so that if the lab sees a pH response, we can tell it came from the histidines.
- **Epitope** — the patch on the target that the binder actually touches.
- **On-site contacts** — how many of those touches land on the specific patch we were aiming for.
- **Pearson correlation (r)** — how well two sets of numbers track each other. 1.0 is perfect, 0 is no relationship.
- **PDB 2TNF** — the published experimental structure of mouse TNF-alpha we built the mouse target from.
- **Liability repair** — small edits that make a sequence easier to manufacture, such as removing a lone cysteine or breaking up a long repeat.

### What is next, and what could go wrong

**Next:** the job running now does two things. It re-examines nine designs the earlier screen skipped because of an arbitrary cutoff — one of them, `tnf2_w0_001`, is in our current submission and looks good on mouse. And it stops the manufacturing-repair step from charging more than the problem is worth. When it finishes I rebuild the 20-design set and we submit. **Deadline is Sunday 18 October, 23:59 AoE — about a day and a half.**

**What is unverified or could be wrong:**

- **No design in this project has been physically tested.** Every number is a prediction from AlphaFold2 against a static structure. The pH switch in particular has never been measured by anyone here — it is a geometric argument plus textbook chemistry about histidine. The lab results will be the first real test.
- **i_ptm is not binding strength.** A design scoring 0.72 is not "twice as good" as one scoring 0.36. It is confidence that the model thinks they stick.
- **The headline correlation finding rests on 11 data points.** The effect is big enough to act on and too small to publish. If anyone asks, that denominator has to come with it.
- **The improvement from the new rule is modest in count.** On a like-for-like comparison the old method yielded 5 usable switches and the new one yields 8. The `tnf2_w0_007` rescue is better evidence that the fix is real than the count is.
- **We optimised for keeping binding, not for switch strength.** The new rule favours single histidines because they disturb less — but a single histidine makes a weaker switch. I have manually overridden this where a two-histidine version also survived, but it means the automatic rule and the competition's top-ranked objective are not perfectly aligned.
- **Two of the ten control twins are not histidine-free** — their parent sequence already contained one. For those two, a pH response in the lab would be harder to attribute cleanly.
- **One design in the proposed set, `tnf2_w0_001`, is carried on older numbers** and was never re-examined by the improved method. The job running now should fix that.
- **I could not check the running job this turn** — the browser connection to Kaggle dropped. Its result is not yet known.
