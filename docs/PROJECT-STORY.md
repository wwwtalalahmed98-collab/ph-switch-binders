# Project story

A running, plain-English record of what this project is doing and why. Newest entry at the top.

| # | Date | Checkpoint |
|---|---|---|
| 7 | 2026-10-24 | The optimiser found the hole in our target before it found the binding site |
| 6 | 2026-10-23 | Depth is not the signal, and we have never once looked at a design |
| 5 | 2026-10-18 | Submitted, and what we would tell the next person |
| 4 | 2026-10-17 | One good histidine out of thirty-five |
| 3 | 2026-10-17 | Submitted, then asked what is actually weak |
| 2 | 2026-10-17 | The set is built and verified; one click left |
| 1 | 2026-10-17 | Measuring the mouse instead of guessing it |

---

## 7. 2026-10-24 - The optimiser found the hole in our target before it found the binding site

### Where we are

The easy half of challenge 3 worked and the hard half failed, for a reason that is our
fault and fixable. We have about a day and a half and one more run's worth of GPU, and I
want to spend it on the design nobody has tried rather than on more of what already works.

### What this project is for

Design small proteins that stick to a chosen human protein, submit at most 20 per
challenge, and have a laboratory test whether they really stick. Challenge 3 is the GLP-1
receptor and it has two requirements, not one: the design must BIND the receptor and it
must ACTIVATE it. A design that binds but does not activate fails the second objective,
and the winner has to meet both.

### Why this step

Two jobs ran overnight on separate GPUs.

The **landing-pad binders** worked: 28 designs, confidence 0.57-0.72, almost all contacting
the groove we aimed at, none of them greasy. These are good molecules.

The **threaders** failed: 1 of 11 is usable. And the way they failed matters more than the
fact that they did.

### What we actually did, and what we found

To make the receptor affordable to design against, we cropped it - kept the part near the
hormone and threw the rest away. Below about 45 angstroms of depth, our cropped target
holds **4 atoms where the real receptor holds 211**. We left a hole.

Ten of the eleven threading designs found that hole and dived into it. In the real
receptor that space is solid protein and those designs would collide with it.

The number worth remembering: across the 11 designs, the model's own confidence correlated
**-0.963** with landing on the residues the natural hormone actually touches. Every design
scoring above 0.66 touched at most 3 of them. The only two that landed correctly scored the
two LOWEST confidences in the set. Ranking by confidence - the obvious thing to do - would
have picked six wrong designs and thrown away the one that works.

This is not the prediction model being unreliable. It did exactly what we asked: find a
sequence that packs confidently against *the target we gave it*. An empty cavity is the
easiest place in any target to pack confidently. We built the flaw; it found it.

### The thing we never tried

Looking at why the natural hormone succeeds: it gets its **grip** from the extracellular
domain and its **activation** from threading the channel. One molecule, two jobs.

We built either/or. Short threaders aimed only at the channel. Folded binders aimed only at
the groove. Nobody required that separation - I introduced it, and it is why our binders
cannot activate and our threaders cannot grip.

The design that matches the mechanism is a longer chain aimed at **both** sites at once.
It has not been attempted, by us or, as far as the published brief suggests, as the
recommended approach.

### The diagram

The rendered diagram in chat shows the receptor drawn three times with the three design
shapes against it: a groove binder that can only bind, a threader that can only reach, and
a two-domain design that does both. The third panel is what the remaining GPU is for.

```mermaid
flowchart TD
    A[GLP-1R: must bind AND activate] --> B[Groove binder<br/>28 made, working]
    A --> C[Threader<br/>1 of 11 usable]
    A --> D[Two-domain design<br/>not yet attempted]
    B --> E[Objective 1 only<br/>this is exendin 9-39]
    C --> F[Failed: 10 of 11 fell<br/>into a hole in our crop]
    F --> G[Seal the target:<br/>192 to 283 residues]
    D --> G
    G --> H[One more run,<br/>both jobs, sealed target]
    style G fill:#993C1D,stroke:#712B13,color:#ffffff
    style H fill:#2563eb,stroke:#1e40af,color:#ffffff
```

"Now" is the blue box.

### What you have now

Working:

- 28 landing-pad binders, measured and filtered. These are the floor: even if everything
  else fails, we have a submittable set.
- [glp1r/results/THREADER_FAILURE.md](glp1r/results/THREADER_FAILURE.md) - the full table,
  the depth-band comparison, and the -0.963 correlation.
- [glp1r/build_submission.py](glp1r/build_submission.py) - the selection rule, written
  before the designs existed, with its one amendment dated and explained.
- Two interactive viewers, in incoming/, for looking at designs rather than reading tables.
- The repository is public and presentable: README with a figure, licence, citation.

Not done:

- The two-domain designs do not exist. Nothing below has been run yet.
- No GLP-1R design has been novelty-checked.
- Nothing tests whether any design ACTIVATES the receptor. Nothing we have can.

### Decisions made, and why

| Chose | Chose against | Why |
|---|---|---|
| Seal the target and re-run | Filter around the hole and accept 1 threader | The hole caused the failure; filtering treats the symptom and leaves objective 2 unaddressed |
| Spend the last GPU on two-domain designs | More landing-pad binders | We already have 28 and need ~10. More of them adds no diversity and cannot activate |
| Validate the sealed target before using it | Launch the long run directly | Twice now a 15-minute check has saved a multi-hour run. The crop change could break the docking it is meant to fix |
| Keep the 28 binders as the floor | Bet everything on the new run | If the two-domain job fails we submit what we would have submitted anyway, minus nothing |
| Re-run threaders on the sealed target | Drop threaders entirely | Same job, two targets, one variable changed - that is a controlled experiment, and the result is worth reporting either way |

### Jargon decoder

- **crop** - cutting a large protein down to the part you need, so the computation is affordable.
- **two-domain mechanism** - the hormone grips one part of the receptor and pushes into another; both are needed.
- **objective 1 / objective 2** - the challenge's two requirements: bind, and activate.
- **confidence (i_ptm)** - the model's own estimate that two molecules fit together. Not a measurement, and this week it pointed the wrong way.
- **correlation** - how reliably two numbers move together. +1 perfectly together, 0 unrelated, -1 perfectly opposite. Ours was -0.963: almost perfectly opposite.
- **on-site** - how many of the receptor residues the real hormone touches are also touched by a design.

### What is next, and what could go wrong

Seal the target, check it still works, then run two jobs overnight: two-domain agonists and
threaders. Select in the morning across three design classes, inspect them in the viewer,
check novelty, write the methods, submit.

What could go wrong, honestly:

- The two-domain designs may simply not fold or dock. It is a harder problem than either
  half, and we will have one attempt at it.
- Sealing the target may disturb the groove docking that currently works. The validation
  step exists to catch that, and if it fails we stop and submit the binders.
- Even a perfect two-domain design is unscreened for activation. We cannot test the thing
  the challenge is actually asking for; we can only build the shape that makes it possible.
- The GPU budget allows one run. There is no second attempt before the deadline.

## 6. 2026-10-23 - Depth is not the signal, and we have never once looked at a design

### Where we are

Challenge 3 opened, we tested our idea before building on it, and the test said the idea
does not work the way we hoped. Separately: we have produced roughly 250 predicted
structures across this project and not once put one on screen.

### What this project is for

Design small proteins that stick to a chosen human protein, submit at most 20 per
challenge, and have a laboratory test whether they really stick. Challenge 2 (TNF-alpha)
is submitted and closed from our side. Challenge 3 is the GLP-1 receptor, closing
Sunday 25 October 23:59 AoE - about two and a half days from now.

### Why this step

GLP-1R is not like the previous targets. It is a receptor that has to be switched ON, and
the switch is mechanical. The hormone's back half grips a landing pad on the outside
(the *extracellular domain*); its first eight residues then thread about 40 angstroms down
a channel into the body of the receptor, and that threading is what fires it. Cut those
eight residues off and you get exendin(9-39) - a real drug that grips the landing pad
perfectly and switches the receptor OFF.

So the obvious design, and the one the challenge page suggests, is a blocker. The thing
that wins is a threader. That meant we needed a way to tell threaders from blockers
before spending the window designing any.

In the TNF challenge we picked our designs using a geometric rule nobody had checked, and
only afterwards measured that the rule pointed +0.24 the wrong way. So this time the check
came first.

### What we actually did

Nature has already run the experiment for us, so we re-ran it in the computer. We fed the
prediction program four peptides against the receptor and measured how deep each one
reached:

| Peptide | What it is in real life | Depth reached | Confidence |
|---|---|---|---|
| GLP-1(7-36) | the real hormone, switches it ON | 34.1 A | 0.73 |
| GLP-1(9-36) | the shape of a drug that switches it OFF | 33.3 A | 0.70 |
| scrambled | same amino acids, shuffled: meaningless | 37.0 A | 0.35 |
| poly-alanine | 30 identical residues: a floor | 30.4 A | 0.18 |

The reference depth measured from the real crystal structure is 33.2 A, so the program did
reproduce the deep binding mode. But it put the blocker in the same place as the hormone -
0.8 A apart, when the measurement noise between two identical re-runs in the TNF work was
up to 0.3 units on a similar score. And the scrambled nonsense sequence reached *deeper
than the real hormone*.

**So depth is not a usable selection rule.** It does not separate on from off, and it does
not even separate sense from nonsense. Confidence does separate sense from nonsense
(0.73 and 0.70, against 0.35 and 0.18), which tells us the model is not broken - it is
simply indifferent to the thing this challenge is about. That is the expected behaviour
of a structure predictor: it has no notion of a receptor being switched on.

This is the second time a geometric proxy of ours has failed a validation test. The
difference is that this time it failed before we spent twenty submission slots on it.

### The second thing: we have no pictures

A separate gap, raised this session. Every structure we have predicted - the 20 submitted
TNF designs, the 35 histidine placements, the GLP-1 controls - exists only inside a Kaggle
notebook's temporary output. Nothing has been downloaded, and nobody has looked at a single
one. We have been selecting designs entirely from numbers in tables.

That is a real weakness, not just a presentation problem. Looking at a structure catches
things a table hides: a binder gripping with a floppy tail instead of a helix, an interface
that is one contact pretending to be ten, a histidine we called "buried" sitting in a gap
full of water.

### The diagram

The rendered diagram in chat shows the three test peptides in the receptor channel, all
landing at the same depth. The point of the picture is that the three bars end level
with each other while their confidence numbers are far apart.

```mermaid
flowchart TD
    A[GLP-1 receptor<br/>two-part mechanism] -->|grip the landing pad| B[Extracellular binder]
    A -->|thread the channel| C[Activating agonist]
    B --> D[Blocks the receptor<br/>this is exendin 9-39]
    C --> E[Switches it on<br/>what actually wins]
    F{Can the predictor<br/>tell B from C?} -->|tested 23 Oct| G[No: 0.8 A apart,<br/>and nonsense scored deeper]
    G --> H[Depth cannot be our<br/>selection rule]
    style G fill:#993C1D,stroke:#712B13,color:#ffffff
    style H fill:#2563eb,stroke:#1e40af,color:#ffffff
```

"Now" is the blue box: we know what we cannot select on, and have not yet decided what
replaces it.

### What you have now

Working:

- [kaggle/cellGLP1_setup.py](kaggle/cellGLP1_setup.py) - builds the receptor target. Ran
  clean: 192 residues kept, all 52 hormone contact points retained, in 1.9 minutes.
- [kaggle/cellGLP2_control.py](kaggle/cellGLP2_control.py) - the four-peptide control test
  above. Open it to see the decision rule written down *before* the numbers came in, which
  is why the PARTIAL verdict is trustworthy rather than reasoning backwards from results.
- [glp1r/scope_6x18.py](glp1r/scope_6x18.py) - measures the hormone's depth profile and
  contact split (30 contacts inside the bundle, 22 on the landing pad) from the real
  structure.

Not done, and not started:

- No GLP-1R designs exist. Zero of 20 slots filled, with two and a half days left.
- No design structure has ever been downloaded or viewed.
- No physics has been applied to anything in this project. Every number we have is either
  a prediction confidence or a geometric measurement.

### Decisions made, and why

| Chose | Chose against | Why |
|---|---|---|
| Validate the depth proxy first | Design 13 agonists and validate later | Exactly the TNF mistake; one 5-minute run avoided repeating it |
| Nature's own control (hormone vs truncated drug) | A synthetic benchmark | The on/off pair exists as real molecules, so the test is about reality, not our assumptions |
| Report PARTIAL as a failure of the proxy | Report "deep mode recovered" as a success | Both are true; the useful half is that it cannot rank candidates |
| Keep the scrambled and poly-alanine floors in the test | Just the two real peptides | The scrambled result is what proves depth is unusable; without it we would have called 0.8 A "a small but real separation" |

### Jargon decoder

- **agonist** - something that switches a receptor on.
- **antagonist / blocker** - something that occupies the site and switches it off.
- **extracellular domain** - the part of the receptor sticking out of the cell; the landing pad.
- **seven-helix bundle** - the part buried in the cell membrane; the channel that must be threaded.
- **exendin(9-39)** - an actual drug: the hormone with its first eight residues removed, which binds and blocks.
- **confidence (i_ptm)** - the predictor's own estimate of whether the two molecules really fit together. Not a measured affinity.
- **angstrom (A)** - a ten-billionth of a metre. The channel is about 40 of them deep.
- **proxy** - a number we can compute standing in for the thing we actually want, which we cannot compute.

### What is next, and what could go wrong

Two decisions are open and they compete for the same two and a half days.

1. **What to submit for Challenge 3.** Depth cannot rank agonist attempts, so any agonist
   we submit is an argued long shot, not a screened candidate. The honest options are a
   majority of landing-pad binders (which our pipeline does well and which will measure
   reliably in the lab) plus a declared minority of threading attempts, or dropping the
   agonist idea entirely.
2. **Whether to spend time on seeing and simulating the designs.** Costs no GPU for the
   viewing part, and it is overdue. The simulation part does cost GPU, and we need the GPU
   for Challenge 3 this week.

What could go wrong: the agonist attempts will almost certainly not activate anything. The
argument for submitting a few anyway is that the selection prompt rewards well-argued long
shots and untested methods - but that is a judgement about the judging, not about biology,
and it should be stated as such.

## 5. 2026-10-18 — Submitted, and what we would tell the next person

### Where we are

The final set is in, it is the one that counts, and the competition is done from our side.

### What was submitted

**Submission 3**, marked DESIGNATED so it supersedes the two before it: twelve acid-OFF switches and eight histidine-free controls, all twenty passing the portal's automatic novelty check. Alongside it went the methods write-up and four results files, including the raw tables behind both negative findings.

The set differs from what we had banked in eight of its twelve switch sequences. Every chosen histidine now faces one of the two anchor points that exist in *both* human and mouse, and **none of them sits in open solvent** — against four of twelve in the previous entry, where the histidine had nothing to push against and could not have worked.

### The honest shape of what we are entering

We are entering a set of designs we believe will mostly fail, with an account of why that is informative anyway.

**The designs are weak.** Median fold confidence 0.693, none above 0.80, where this field normally throws away anything below 0.85 before ordering. Every one carries a sequence flaw of the kind that predicts poor manufacturing. We skipped the standard clean-up step that would have fixed this, because doing it properly would have thrown away every measurement we had made, with a day left.

**The switch is still not proven.** Across all 35 ways we knew how to build it, one placement showed a clear effect and 34 did not. The chain from what we measured to an actual acid-release has two unproven links in it. We said so in the submission rather than letting the numbers imply more than they support.

**My honest estimate of the odds** has not changed: a few percent that any single design binds, and if it binds, perhaps one in four that the switch is clean. Across twelve, maybe a one-in-five chance of one real hit.

### What we would actually tell someone doing this next

Two findings, both negative, both the same shape — a cheap proxy we were trusting turned out to carry almost no information:

1. **Sequence similarity does not predict cross-species binding.** We ranked mouse compatibility by how similar the contact points looked between species. Correlation with the real thing: +0.272. The design it ranked first had the worst actual mouse binding in the set.

2. **Distance does not predict the acid switch.** We placed histidines by how close they sat to a positive spot on the target. Correlation with whether that spot actually rejects a positive charge: +0.24. Whether *gripping buries* the histidine: −0.58. We had been selecting on the wrong variable the whole time.

The second one only came to light because we went looking for evidence that our own method worked, and did not find it.

### The diagram

```mermaid
flowchart TD
    A[89 designs generated] -->|31 strong on human| B[Tested on mouse directly]
    B -->|20 of 31 bind both| C[35 switch placements<br/>that hold both species]
    C -->|does the spot reject<br/>a positive charge?| D[1 of 35 clears the noise]
    C -->|does gripping bury<br/>the histidine?| E[Rank on this instead]
    E --> F[12 switches + 8 controls<br/>SUBMITTED, designated]
    D -.->|informs, does not decide| F
    style F fill:#15803d,stroke:#14532d,color:#ffffff
    style D fill:#b45309,stroke:#78350f,color:#ffffff
```

Green is in and counting. Amber is the uncomfortable number we reported rather than hid.

### Decisions made, and why

| We chose | Instead of | Because |
|---|---|---|
| Replacing the banked set | Leaving it alone | 10 of 12 placements were beaten by another position in the same design, and 4 had a histidine that could not work |
| Twelve single histidines | Six pairs | Every pair had a dead partner, so the two-histidine chemistry argument did not apply to any of them |
| Passing over two higher-burial placements | Taking the top of the ranking | They faced a cation with no mouse equivalent, and mouse ranks above affinity |
| Leading the write-up with two negative results | Leading with the designs | They are the most transferable thing we produced, and method novelty is a selection criterion |
| Stopping measuring | Spending the last day on more runs | The remaining questions needed work that would invalidate what we had, not add to it |

### What is next, and what could go wrong

**Nothing from us.** Submissions close 19 October 11:59 UTC. Selection happens afterwards, by a Claude-based workflow reading the methods alongside the designs, and any chosen designs go to Adaptyv's lab. Every measurement, including the failures, gets published.

**What could still go wrong:** all of it, in the ordinary way. The designs may not express, may not bind, and the switch may not fire. We have said as much in the submission. The part we would stand behind regardless is the two correlations — those are measurements, and they would hold whether or not a single one of our designs works.

---

## 4. 2026-10-17 — One good histidine out of thirty-five

### Where we are

We measured all 35 ways we know how to build the switch, found that only one of them has any real support, and now have to decide whether to replace what is already entered.

### What this project is for

Anthropic's protein design competition, Challenge 2, Track 3. A small protein that grabs human TNF-alpha at normal body acidity, lets go when slightly acidic, and does the same on the mouse version. Ranked on the acidity switch first, mouse second, binding strength third.

### Why this step

We had submitted twenty designs chosen on a rule that nobody had ever tested: put a histidine within a certain distance of a positively charged spot on the target. The previous run showed that rule barely predicts anything, and that a better one — how much *gripping* buries the histidine — predicts three times better. This step applied the better rule to every placement we have, not just the twelve we happened to send.

### What we actually did

Measured all 35 placements across 14 design families: how buried each histidine gets when the binder grips, whether that spot actually objects to a positive charge, whether the binder folds on its own, and how far each histidine really is from a positive spot once folded.

It ran in **40 minutes**, against 79 minutes for the previous smaller job, because I stopped it rebuilding the prediction model between every single prediction.

### The three things it found

**One placement in thirty-five has real support.** We test each histidine spot by asking what happens if we put a *neutral* amino acid there versus a *positively charged* one. If the spot genuinely rejects a positive charge, the switch has a reason to fire. Across 35 placements, exactly one shows a clear effect: `tnf2_w0_001` Y37H. It is also the only one that reproduced — an earlier run measured the same spot independently and got the same answer. For the other 34, the effect is smaller than the wobble we get re-running the identical calculation, and in **8 of them a positive charge actually helps**.

**Every two-histidine design has a dead partner.** Our argument for using two histidines was simple chemistry: one histidine leaves half the molecules still gripping in acid, two leaves a quarter. But that assumes both are doing work. Ranked on the *weaker* of their two histidines, every pair we have sits near the bottom: the best pair scores 0.100 where the best single scores 0.404, and our top submitted design scores **0.000** — one of its two histidines sits in open water and cannot contribute anything. Across all 35, **seven have a histidine that is doing nothing**.

**I was wrong about the folding.** Last checkpoint I said the designs might only hold their shape when the target is there to hold them, inferred from Proteinbase reporting lower confidence than we did. Measured directly: they score *higher* on their own than when gripping. They are fine. The inference was wrong, and the reason it cost nothing is that I flagged it as a guess and measured it instead of acting on it.

### The bad news

**Only two of our twelve submitted placements are the best their family offers.** The other ten are beaten by a different histidine position in the same design.

**The evidence for the switch is thin and I do not want to oversell the fix.** Burial is a *correlate* of charge rejection, not a measurement of switching, and charge rejection is itself a stand-in for chemistry the prediction model does not actually simulate. There are two unproven links in that chain. What we can say is that we looked for support at 35 placements and found it at one.

**Eleven of 41 histidines face a positive spot that is not one of the two conserved ones** — meaning even if those switches fire in human, there is no guarantee of the same spot existing in mouse.

**One predicted structure is simply wrong**: a histidine sits 0.86 angstroms from a positive charge, which is closer than two atoms can physically be.

### How it works

The hand-and-ball picture again: the binder is a hand gripping a ball, the switch is a spring-loaded spike that pops out in acid and pushes the hand off.

What we learned is that a spike only works if **closing the hand traps it**. A spike sticking out into open air pops out freely and the hand keeps gripping — it costs nothing and does nothing. That is what "zero burial" means, and seven of our spikes are like that, including one of the two on our best design. Two spikes are only better than one if both are trapped.

### The diagram

```mermaid
flowchart TD
    A[35 ways to place the switch] -->|does the spot reject<br/>a positive charge?| B{Real effect?}
    B -->|1 placement| C[Y37H, and it reproduced]
    B -->|34 placements| D[Smaller than our own noise]
    A -->|how much does gripping<br/>trap the histidine?| E{Trapped?}
    E -->|28 placements| F[Something to work with]
    E -->|7 placements| G[Sits in open water,<br/>cannot contribute]
    F --> H[Best single beats<br/>every pair]
    G --> H
    style C fill:#15803d,stroke:#14532d,color:#ffffff
    style G fill:#b45309,stroke:#78350f,color:#ffffff
```

Green is the one placement with evidence behind it. Amber is the seven that cannot work regardless of chemistry.

### What you have now

- [tnf/results/BURIAL_ALL_RESULTS.md](../tnf/results/BURIAL_ALL_RESULTS.md) — all 35 placements with burial, charge, and binding, and the best option for each family. Open it for the table showing which of our submitted choices were beaten.
- [tnf/results/SWITCH_EVIDENCE_RESULTS.md](../tnf/results/SWITCH_EVIDENCE_RESULTS.md) — the run before it, which found the better rule.
- **Submitted and counting:** 12 switches and 8 controls, banked and designated. Unchanged by any of this.

### Decisions made, and why

| We chose | Instead of | Because |
|---|---|---|
| Measuring all 35 placements | Only the 12 we sent | 10 of the 12 turned out not to be the best their family offers |
| Reporting one positive out of 35 | Reporting the best few | Only one cleared the noise floor, and saying more would be inventing signal |
| Measuring the folding claim | Repeating last checkpoint's inference | It was a guess, and it was wrong |
| Running measurements in grouped passes | The previous per-design order | The old order rebuilt the model 152 times for 212 predictions |

### Jargon decoder

- **Burial on binding** — how much of the histidine gets covered up when the binder grips the target. A histidine that stays exposed cannot drive release.
- **Charge test** — predicting the same spot as a neutral amino acid and as a positively charged one. The gap says whether the spot objects to charge.
- **Noise floor** — how much a number moves when you re-run the identical calculation. Ours is about 0.02 to 0.03, so smaller differences mean nothing.
- **Conserved anchor** — a positively charged spot present in both the human and mouse target. A switch built against a human-only spot cannot work in mouse.
- **Henderson-Hasselbalch** — the textbook sum for how much of a histidine is charged at a given acidity. The basis for preferring two histidines over one.

### What is next, and what could go wrong

**Next:** decide whether to replace the submitted set. The window opens 18 October 08:14 UTC and the deadline is 19 October 11:59 UTC, so there is about a day of room.

**The risk either way.** Replacing means trading six two-histidine designs — backed by textbook chemistry — for single-histidine ones chosen on our own unvalidated correlate. Not replacing means knowingly leaving in designs whose histidines provably cannot contribute. The strongest part of the case for replacing does not depend on the shaky chain at all: a histidine with zero burial cannot help, which is geometry rather than prediction, and the replacements are also better on *measured* binding in six of the families.

---

## 3. 2026-10-17 — Submitted, then asked what is actually weak

### Where we are

The designs are submitted and counting, and we have started measuring the one thing we had never measured: whether the acidity switch works at all.

### What this project is for

Anthropic's protein design competition, Challenge 2, Track 3. A small protein that grabs human TNF-alpha at normal body acidity, lets go when slightly acidic, and does the same on the mouse version. Up to 20 designs, ranked on switch first, mouse second, binding strength third.

### Why this step

Last checkpoint the set was built but blocked at the upload. This step banked it, then asked the harder question: not "is the set internally consistent" but "is this good science, and will any of it work in the lab".

### What we actually did

**Changed the set before sending it.** We had 10 switches plus 10 histidine-free twins. A twin exists so that if the lab sees an acidity response we can prove our histidine caused it — but a twin can never win the top-ranked objective itself. Four of our twins carried a histidine of their own, so they could not serve as a clean baseline anyway. We gave those four slots to more switches. Result: **12 switches and 8 controls**, and the number of compromised controls went from 2 to **0**. More shots at the thing that matters most, and cleaner evidence.

**Submitted.** All 20 passed the portal's automatic novelty check. The entry is marked DESIGNATED, so it is the one that counts and it replaces the earlier set in which six of ten switches failed the mouse objective.

**Then looked hard at quality.** Honest summary: the designs are weak as folds. Median pLDDT 0.693, and **none of the 20 reaches 0.80**. Published work in this field filters at 0.85 before ordering anything. Every single design also carries a sequence flaw — long runs of one amino acid, big greasy stretches, repetitive patches — which are the fingerprints of the generator we used and a plausible cause of the protein simply not being manufacturable.

**Found that one of those flaws was ours.** Our cleanup step breaks up runs of repeated amino acids. On one design it turned a run of Q into an alternating QNQNQN pattern. That is not a fix; it is the same problem wearing a different hat.

**The portal handed us a free experiment.** Proteinbase runs its own structure prediction, and reported confidence of 55 to 63 where ours said 60 to 77 for the same designs. The most likely explanation is that it predicted the binder **on its own** while we predicted it **stuck to the target**. If so, our designs may only hold their shape when TNF-alpha is there to hold them — which is a poor sign for making them at all. We did not take this on trust; the run now going measures it directly.

### The bad news

**We have never tested the switch. Not once.** Everything measured so far is binding at normal acidity against two species. No calculation in this project touches acidity. The top-ranked objective has been supported entirely by "we put a histidine near an arginine".

**And the rule we used for that was the wrong rule.** We accepted any histidine between 3.0 and 8.8 angstroms from a positive charge. At the far end of that range, in salt water, the repulsion is a fraction of the energy needed to matter. The quantity that actually predicts the switch is how much *binding buries* the histidine — because in the free binder it is exposed and picks up its charge happily at pH 6, and it is burying that charge against the target that costs the binding. We filtered on the wrong variable, and we are now computing the right one.

**Our whole ranking has no error bars.** Every number is a single prediction. The spread across model versions is normally about the same size as the margins our ranking depends on, so some of our ordering is noise and we cannot yet say which part.

**The biggest known fix is one we chose not to attempt.** A redesign pass (ProteinMPNN) is the standard step that raises fold quality. We tried it earlier in the competition, used it wrongly, saw it wreck the interfaces, and dropped it. Doing it properly now would regenerate every sequence and throw away all 35 switch measurements, with a day left. We named it in the submitted methods as the largest gap rather than hiding it.

**Honest odds.** My estimate, and it is only an estimate: a few percent that any one design binds at all; if it binds, maybe one in four that the switch is clean. Across twelve switches, perhaps a one-in-five to one-in-four chance of a single real hit. That is roughly what the organisers expect — they screen hundreds knowing most fail — but it means our realistic route to mattering is the method and the honest negative result, not raw binding strength.

### How it works

Keep the hand-and-ball picture: the binder is a hand gripping a ball, the switch is a spring-loaded spike on one finger, flat normally and popped out in acid.

The correction this time is about *where* the spike has to sit. We had been measuring how close the spike is to a bump on the ball. What actually matters is whether closing the hand **traps** the spike with nowhere to go. A spike sticking out into open air can pop out freely and the hand keeps gripping. A spike that the grip itself presses into a tight pocket cannot pop out without forcing the hand open. That trapping is what we should have been measuring, and now are.

### The diagram

```mermaid
flowchart TD
    A[10 switches + 10 twins] -->|4 twins had their own<br/>histidine, so could not<br/>serve as a baseline| B[12 switches + 8 clean twins]
    B -->|all 20 pass<br/>novelty check| C[Submitted and DESIGNATED]
    C --> D{But is any of it good?}
    D -->|folds are weak:<br/>0 of 20 above 0.80| E[Named in the methods]
    D -->|switch never tested| F[Run now measuring it]
    F -->|histidine as neutral<br/>vs as charged| G[Does the spot reject<br/>a positive charge?]
    F -->|buried by binding?| H[The right variable<br/>at last]
    style F fill:#2563eb,stroke:#1e40af,color:#ffffff
    style C fill:#15803d,stroke:#14532d,color:#ffffff
```

"Submitted" (green) is banked and counting. "Run now measuring it" (blue) is where we are.

### What you have now

**Things that work and have produced results:**

- [tnf/results/SWITCH5_RESULTS.md](../tnf/results/SWITCH5_RESULTS.md) — the run behind the submitted set, and what it replaced.
- [submission/METHODS_TNF_V3.md](../submission/METHODS_TNF_V3.md) — what reviewers read. Open it to see the limitations stated plainly, including the ones found this turn.
- [tnf/build_submission_v3.py](../tnf/build_submission_v3.py) — rebuilds the exact 20 from the run output and refuses any switch that differs from its twin anywhere but the histidines.

**Submitted and counting:** 12 switches and 8 controls on Proteinbase, marked DESIGNATED. The design files stay out of the public repository until the competition closes.

**Running now:** [kaggle/cellTNF20_switch_evidence.py](../kaggle/cellTNF20_switch_evidence.py), measuring four things — whether each histidine position rejects a positive charge, error bars over five model versions, whether each binder folds on its own, and how much binding buries the histidine.

**Not done:** the redesign pass that would fix fold quality.

### Decisions made, and why

| We chose | Instead of | Because |
|---|---|---|
| 12 switches + 8 controls | 10 + 10 | A twin cannot win objective 1, and four of ours carried a histidine so could not even serve as a baseline |
| Submitting before running the new tests | Waiting for better numbers | The banked set had six of ten failing mouse. Banking an improvement costs nothing, and the portal allows another submission in 24 hours |
| Testing the histidine as neutral vs as charged | Trusting the geometric argument | It is the only falsifiable test of objective 1 available, and it can tell us a design's switch claim is unsupported |
| Measuring burial on binding | Keeping the 3.0-8.8 angstrom distance rule | At the far end of that window the repulsion is too small to matter; burial is the quantity that actually drives the switch |
| Not attempting the redesign pass | Fixing fold quality properly | It regenerates every sequence and invalidates all 35 switch measurements, with a day left. Named as the largest gap instead |
| Stating the weaknesses in the submitted methods | Presenting the set at its best | A Claude-based workflow reads these when choosing what to make, and an undisclosed flaw a reviewer can spot is worse than a disclosed one |

### Jargon decoder

- **pLDDT** — the prediction model's confidence in the shape it drew, 0 to 1. Below about 0.8 for a designed miniprotein is a warning sign.
- **i_ptm** — confidence that two proteins stick together. A confidence score, not a measured strength.
- **ProteinMPNN** — a tool that rewrites a protein's sequence to fit a shape better. The standard quality step we skipped.
- **Novelty check** — the portal's automatic test that a design is genuinely new rather than a copy of a known protein. Needs 3 of 4; all 20 passed.
- **DESIGNATED** — the portal's label for which of several submissions actually counts.
- **Debye screening** — salt in water cancels out electrical attraction and repulsion over distance. Why a charge 8.8 angstroms away does almost nothing.
- **Avidity** — grabbing something in several places at once is much stronger than the sum of the grips. Works against us here, because partial loss of grip in acid may still read as "still holding".

### What is next, and what could go wrong

**Next:** read the run, then decide whether it changes the set. The portal allows one more submission in about 24 hours, and the deadline is Sunday 18 October 23:59 AoE, so there is comfortable room for one upgrade.

**What could go wrong:**

- **The charge test could come back negative across the board.** If histidine positions turn out to tolerate a positive charge, our switch claim is unsupported and we would have to say so. That is the point of running it, but it is a real possible outcome.
- **The test is weak.** The prediction model is not a reliable energy calculator, so only large effects mean anything. A small shift tells us nothing either way.
- **The monomer result could confirm the bad reading** — that our binders do not fold without the target. We would then know our designs are unlikely to express, with no time to fix it.
- **Everything from the previous checkpoints still stands**, led by the fact that the acidity switch is reasoned rather than computed, and that nothing has been physically tested.

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
