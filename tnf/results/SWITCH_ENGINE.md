# The acid-OFF switch engine, and two corrections it exposed

Challenge 2 inverts Challenge 1: bind at pH 7.4, release at pH 6.0. Histidine now
faces a target **Arg or Lys** instead of a carboxylate. Neutral at 7.4 it is
tolerated; protonated at 6.0 it collides with a charge already present and the
interface gives way.

## Calibration: measure what nature avoids

Challenge 1 asked where nature *puts* a histidine next to a carboxylate. The acid-OFF
question is the mirror image — where does nature *avoid* putting one next to a cation?
Those are precisely the places where protonation costs something.

Measured over 17 structures, 278 His–carboxylate and 279 His–cation closest approaches:

| closest sidechain N → partner | p5 | median | within 4.0 Å |
|---|---|---|---|
| His – carboxylate | 2.81 Å | 5.93 Å | **22.3%** |
| His – cation | 3.57 Å | 7.74 Å | **9.7%** |

Nature places a histidine within 4 Å of an acid **2.3× more often** than within 4 Å of
a cation. That ratio is the signal the switch exploits.

A second result was not expected: the **Cβ distance distributions of the two cases are
nearly identical** (median 6.54 vs 6.39 Å). Contact geometry does not distinguish
attraction from repulsion — only frequency does. So the distance window transfers
between the two directions and only the partner and the scoring shape change. Better
to state that than re-derive a window and pretend it was new.

The scoring shapes differ because the physics does. Repulsion rises monotonically as
charges approach, so acid-OFF quality increases as the contact tightens, weighted by
burial since solvent screens charge.

## Validation on a switch the engine was never shown

**4N0U** is human FcRn bound to an IgG Fc — the textbook pH-dependent interaction. IgG
binds FcRn in the acidified endosome and releases at the neutral pH of plasma, through
Fc **His310** and **His435**. It is an acid-ON switch in a system unrelated to anything
used for calibration.

Scored on geometry alone (the mutation-preference term is meaningless for a residue
that is already histidine, and would penalise the thing under test):

```
acid_on    11 candidate contacts across 10 positions
           His310  rank 3 of 10      H310 -> E115(A)  Cb 3.09 A
           His435  rank 9 of 10
acid_off    1 candidate contact
           His310  not found         His435  not found

acid-ON  identifies 2 of the 2 known switch residues
acid-OFF identifies 0 of them  (correct - this complex is acid-ON)
```

Both known residues recovered, in the right direction, with the opposite mode
correctly silent.

## Two corrections this exposed in the Challenge 1 engine

**1. The distance window was too narrow at both ends.** Challenge 1 used Cβ 3.5–7.0 Å.
The control has His310 at **3.09 Å** and His435 at **7.46 Å** — the textbook switch
would have been missed at *both* residues. The calibration's own percentiles give p1
≈ 3.0 and p95 8.70; the window is now **3.0–8.8 Å**. The original figures were set by
eye from the measured contacts rather than from their percentiles.

**2. The quality function was backwards for attraction.** Challenge 1 scored acid-ON
as a symmetric peak at the median contact distance, which treats a *tighter* salt
bridge as worse than an average one. What limits an attractive pair is whether the
imidazole can reach, so quality is now flat inside comfortable reach (≤6.5 Å) and
falls only as reach runs out. Orientation was already handled separately.

Both corrections affect Challenge 1's results, not just Challenge 2's. The over-strict
window is a plausible contributor to how few switch sites that search found — it
reported zero across 14 designs before the carboxylate set was widened, and the window
was narrowing the field further at the same time.
