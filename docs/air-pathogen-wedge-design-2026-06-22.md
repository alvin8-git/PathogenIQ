# Design: the minimum wedge for air-sample pathogen detection

**Date:** 2026-06-22 · **Mode:** intrapreneurship/research · **Status:** approved (office-hours)

## Problem

We want to detect **possible pathogens** in low-biomass air samples and explicitly
**not** catalog environmental microbes. The pipeline has been accreting tools
(PhiX, VFDB, spike-in quant, and a heavy assembly/MAG arm needing ~100 GB of
GTDB-Tk/CheckM DBs). The concern: scope explosion. The question: what is actually
necessary for the air-pathogen wedge?

## Core finding

The biggest exploder — the **assembly/MAG arm (MEGAHIT → MetaBAT2 → CheckM →
GTDB-Tk)** — exists to recover **novel/environmental organisms**, which is the
*opposite* of the stated goal. Naming Sphingomonas/Methylobacterium is
environmental ecology (the Jeilu et al. survey), not pathogen detection. The
heaviest dependency in the pipeline serves the goal we said we don't want.

**The minimum wedge for air pathogen detection is small, and we already have most
of it:**

| Core piece | Status | Why it's in the wedge |
|------------|--------|-----------------------|
| Targeted pathogen DB (tier-1, ~110 + biothreat) | have | Narrow on purpose — "pathogens not environmental microbes" |
| **NTC background subtraction** (`background.py`) | have | THE load-bearing piece for air; the kitome out-masses signal |
| Evidence grading + NTC tier cap (`report.py`) | have | Turns hits into A/B/C/X with controlled confidence |
| Breadth-of-coverage gate (`coverage.py`) | have (needs wiring) | Real genome-wide signal vs one-locus contamination artifact |
| **Novelty flag** (NEW, `novelty.py`) | to build | Stay un-blind to a novel pathogen without the 100 GB arm |

Everything else (assembly, GTDB taxonomy, diversity/ordination, HUMAnN) is
**environmental-survey tooling** — a separate optional module, not the core.

## Decision

1. **Scope = known pathogens + a cheap novelty flag** (not full novel-organism
   discovery in the core).
2. **Novelty flag = Approach A:** run the existing **Kraken2 Standard-8 DB (8 GB)**
   over the non-host reads, compute the **unclassified fraction** (reads matching
   nothing in RefSeq = candidate novelty); flag when it is high *and* complex.
   Reuses `benchmark.py`'s Kraken parser. The 8 GB DB is the price of not being
   blind — vs the 100 GB GTDB-Tk arm.
3. **Assembly/MAG arm → Tier-2, triggered.** It stays CLI-gated (`--assemble`),
   demoted from the always-on core. It runs only when the novelty flag fires (v1:
   manual; auto-escalation is a v2, Approach B).

## Plan 6 re-triage against the wedge

| Item | Verdict |
|------|---------|
| #1 PhiX removal | **Core** (cheap QC correctness) — done |
| #6 VFDB virulence | **Core, and important** — separates a pathogen from a same-genus commensal — done |
| #2 spike-in absolute quant | **Core-adjacent** — "is the pathogen above an actionable load?" — done |
| #3 assembly/MAG | **Demote to Tier-2 triggered** (built, CLI-gated; not core) |
| #4 BHI enrichment | Wet-lab, optional — defer |
| #5 strain ID | Defer (not core pathogen detect) |
| #7 diversity/ordination | **Out of core** — environmental ecology |
| #8 GTDB taxonomy | Only inside the triggered assembly arm |
| #9 overrepresented-seq trim | Minor QC — optional |

## New work for the wedge (small)

1. **`novelty.py`** — Kraken2 Std-8 over non-host reads → unclassified fraction +
   complexity guard → `NoveltyFlag` (fraction, n_unclassified, triggered:bool).
   Surface in the report; calibrate the threshold on the air NTCs. CLI: `--novelty`.
2. **AIR specimen type** — the pipeline requires a `SpecimenType`; air has none
   (the concordance run used `blood` as a stand-in). Add `AIR` with air-appropriate
   read thresholds and air contaminant priors (Cutibacterium/Sphingomonas/
   Methylobacterium/Ralstonia kitome). This is a real gap, not a nice-to-have.
3. **Air NTC as the standard background** — already supported via `--background`;
   document the air-NTC build (`scripts/11` on `air-aircraft-ntc`) as the
   recommended air workflow. NTC is mandatory for air.

## Why this is the right wedge (premises, agreed)

- The wedge is **interpretation, not detection** — detection is commodity.
- A **narrow targeted DB is a feature** — you want to ignore environmental reads.
- **NTC is mandatory for air** — without a control, low-biomass air is
  uninterpretable.
- The novelty flag's broad reference is the **8 GB Kraken2 DB we already have**,
  not GTDB-Tk's 100 GB.
- Assembly/GTDB is **Tier-2 discovery, triggered**, not core.

## Alternatives considered

- **B (formal two-tier, auto-escalation):** cleaner long-term but more wiring and
  re-introduces the gated GTDB dep — premature before the flag proves useful. The
  natural v2.
- **C (DB-free novelty from unmapped-pool stats):** lightest footprint but cruder
  and calibration-risky; could miss a coherent novel organism a classifier catches.

## The assignment (do this next, real-world)

**Prove the wedge on real air before building the flag.** The air concordance run
(in progress, `bqm3jgyqh`) runs targeted detection + air NTC + grading on the 6
aircraft filters. When it lands, answer two questions from the JSON:
1. Of the paper's **in-DB** taxa (*E. coli*, *P. aeruginosa*), how many are
   detected and **gradeable** after air-NTC subtraction?
2. **How much kitome did the air NTC remove** (taxa dropped as background)?

If NTC + grading cleans air down to a short, interpretable pathogen list, the
wedge is proven and `novelty.py` is the only new core code. If it doesn't, the
NTC model (not the novelty flag) is where the next work goes.

## Out of scope (separate module, do not build into core)

De novo environmental cataloging, MAG taxonomy at scale, functional/diversity
profiling, strain phylogenetics. These are a *survey* product, not the
*pathogen-detection* wedge.
