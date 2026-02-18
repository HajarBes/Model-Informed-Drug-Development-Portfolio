# AUCR Definition and Interpretation

## Metric 1: Stratum-Level Median Shift

**Definition:**

    AUCR_shift = median(AUC at weight W) / median(AUC at 70 kg)

**Denominator:** Median AUC across N=1000 simulated subjects at 70 kg
reference weight (7.5 mg dose, allometric scaling applied).

**What it measures:** The shift in typical (population-median) exposure
attributable to body weight. This isolates the weight effect by comparing
stratum-level summary statistics. It is a single number per weight
stratum and approximately equals the theoretical allometric prediction:
(70/W)^0.75.

**Output file:** `aucr_stratum_shift_summary.csv`

## Metric 2: Individual-Level AUCR Distribution

**Definition:**

    AUCR_i = AUC_i / median(AUC at 70 kg)

For each simulated subject i at weight W, the individual AUC is divided
by the population median at 70 kg.

**What it measures:** The combined effect of weight AND between-subject
variability (BSV) on individual exposure relative to a typical 70 kg
subject. This metric reflects what a clinician would observe in practice:
a specific patient may deviate from reference due to both their weight
and their individual metabolic capacity.

**Why 45% of the 70 kg reference group falls outside 0.80-1.25:**
Even at the reference weight, individual AUC varies because of BSV in
CL (~30% CV). The lognormal distribution of CL produces a spread in
individual AUCs where roughly half of subjects fall outside the
0.80-1.25 band around the median. This is not a weight effect — it
is pure pharmacokinetic variability.

**Output file:** `aucr_weight_summary_verified.csv`

## How to Read the Two Metrics Together

| Metric | Tells You | Example at 50 kg |
|--------|-----------|------------------|
| Stratum shift (1.33) | Weight alone increases typical AUC by 33% | Modest, predictable effect |
| % outside band (62.8%) | Most individuals deviate >20% from 70 kg median | Driven by BSV, not weight |

The stratum shift is the **weight effect**. The % outside is
**variability + weight combined**. To assess whether dose adjustment
is needed, focus on the stratum shift. Use the individual-level
distribution to understand what range of exposures a prescriber
should expect.

## The 0.80-1.25 Band

The 0.80-1.25 interval is used here as a **clinical relevance heuristic**
— a rule of thumb for judging whether covariate-driven exposure shifts
are likely to matter in practice. It is borrowed from the bioequivalence
literature (where it has formal statistical meaning for AUC ratios of
geometric means) but applied informally to median exposure ratios across
weight strata.

This is NOT:
- A formal bioequivalence assessment
- A regulatory acceptance criterion
- A statistical hypothesis test

It is a convenient benchmark for clinical interpretation, analogous to
the 80-125% convention used in forest plots across the pharmacometrics
literature.
