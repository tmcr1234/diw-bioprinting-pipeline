# SOP — The printability index: one definition, and what it can and cannot decide

**Version** v1 · 2026-08-26 · Status: `active`
**Scope** Stage 4 (DIW print validation). Applies to every use of the Ouyang
printability index in this lab, in code and in manuscripts.

This SOP exists because the same index was defined two incompatible ways in
two sections of one manuscript, and was separately used to make a claim its
replicate count could not support. Both were caught by external audit, late.
Neither is a hard problem; both are a documentation problem.

---

## 1. The canonical definition — there is exactly one

```
Pr = L² / (16·A)
```

where, **per inter-strand pore**, `L` is the pore perimeter and `A` is the
pore area. A perfectly square pore gives `Pr = 1`; under-gelation rounds the
pore and drives `Pr` below 1; over-gelation squares the corners and drives it
above 1.

**The implementation is the definition.** It lives in exactly one place:

> `Export/04_Printability/Printability/metrics_pore.py`

with a unit test in `metrics_pore`'s test module fixing a known pore at
`Pr ≈ 1.06`. A manuscript must **cross-reference** that file, never restate
the formula from memory in a Methods section and again in Results. The two
restatements in the incident were `expected/measured` in one section and
`measured/expected` in the other — reciprocals of each other. The wrong one
would have flipped the sign of the die-swell conclusion that rested on it.

**Rule.** If you find yourself typing `Pr =` into a manuscript, stop and
write "computed per Eq. (n) using the implementation in `metrics_pore.py`"
instead. If a paper genuinely needs the formula inline, it appears **once**,
in Methods, and every later mention refers back to it.

### The ideal value is not 1 for a rectangular target

For a target pore of width `W` and height `H`,

```
Pr_ideal = (2·(W + H))² / (16·W·H)
```

which is 1 only when `W = H`. Reporting `Pr` against 1 for a deliberately
rectangular lattice is a category error. `metrics_pore.py` documents this;
manuscripts must not silently assume the square case.

---

## 2. What the index can decide, and at what replicate count

### The reader's objection, stated as a heading: *"you called one formulation better — can your data tell them apart?"*

In the incident, an `n = 2` per cell design gave a formulation-effect ANOVA
of **p = 0.20**. The index simply could not discriminate the best from the
worst formulation at that replicate count. The manuscript had to be
rewritten across ten versions (v25 → v35) to reframe the index as a
**go/no-go screen** rather than a ranking validator.

That reframe was correct. It was also entirely avoidable, and it exists only
in that manuscript's prose — not anywhere a future student would find it.

### The arithmetic, spelled out

Two formulations, `n` prints each, comparing mean `Pr`. To detect a
difference `Δ` between formulation means with a two-sample t-test at
α = 0.05 and 80 % power, you need roughly

```
n  ≈  16 · σ² / Δ²          per group
```

(the standard `2·(z_{α/2}+z_β)² ≈ 15.7` rounded up). So with a between-print
standard deviation of `σ = 0.05` in `Pr`:

| difference you want to detect, Δ | 16·σ²/Δ² | n per group |
|---|---|---|
| 0.02 (2 % of ideal) | 16·0.0025/0.0004 | **100** |
| 0.05 | 16·0.0025/0.0025 | **16** |
| 0.10 | 16·0.0025/0.01 | **4** |
| 0.20 | 16·0.0025/0.04 | **1** (→ use 3, the floor) |

Read the table the honest way round: **at n = 2 you can only detect
differences of about Δ ≥ 0.20 in `Pr`** — a gross difference, the kind you
can see by eye. Anything finer is below the resolution of the design.

> **Measure your own σ first.** The 0.05 above is illustrative. Compute the
> between-print standard deviation from your own replicate prints of a
> single formulation, then read the table with that number. `aggregator.py`
> already emits `Pr_consensus_std` per print — use it.

### The rule

| replicate count | what you may claim |
|---|---|
| n = 1 | Nothing quantitative. A picture. |
| n = 2–3 | **Go / no-go screen only.** "Formulation X is printable; Y is not." Never a ranking, never "X is better than Y". |
| n ≥ 6 | Ranking claims become defensible for differences of Δ ≈ 0.08 at σ = 0.05. State the power calculation. |
| n from the table | Whatever Δ the table supports at your measured σ. |

### When the replicate count is what it is

Sometimes n = 2 is all the material allows. That is fine, and it does not
invalidate the work — it changes what the work is allowed to say. Then:

1. Report the ANOVA **with its p-value**, and say plainly that it does not
   reach significance.
2. Frame the index as a screen in Methods, before any results are shown —
   not as a retreat in Discussion after a reviewer objects.
3. Report the observed Δ and the Δ your design could have detected, so the
   reader can see the gap themselves.
4. Do not rank. Do not write "best", "optimal", or "superior".

### Before and after, on the incident

| | before (v25) | after (v35) |
|---|---|---|
| role of the index | ranking validator — "formulation C is the best" | go/no-go screen — "C15-SF5.5 is printable" |
| n per cell | 2 | 2 (unchanged — the *data* was never the fix) |
| ANOVA reported | no | yes, p = 0.20, stated as non-significant |
| claim strength | unsupported | supported |
| rounds of revision | — | 10 |

The data never changed. Only the claim did.

---

## 3. Checklist before a printability number leaves this lab

- [ ] `Pr` came from `metrics_pore.py`, not from a formula retyped into a script.
- [ ] The manuscript states the formula **once** and cross-references it thereafter.
- [ ] `Pr_ideal` accounts for a non-square target pore, if the target is non-square.
- [ ] Between-print σ measured from replicates of one formulation.
- [ ] Replicate count checked against the table in §2 for the Δ actually claimed.
- [ ] If n ≤ 3, the text frames the index as a screen and makes no ranking claim.
- [ ] ANOVA p-value reported whether or not it is significant.

---

*Related:* `Export/04_Printability/Printability/metrics_pore.py` (definition),
`aggregator.py` (`Pr_consensus_mean` / `Pr_consensus_std`),
`Export/CLAUDE.md` §4 (slicer convention, which feeds the width the index sees).
