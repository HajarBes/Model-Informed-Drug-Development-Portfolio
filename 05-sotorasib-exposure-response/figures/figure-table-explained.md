
---
  Figures

  Figure 1: dose_exposure_boxplot.png

  What it shows: Side-by-side boxplots of simulated steady-state AUC (left) and Cmax (right) across all five dose
  levels (180–960 mg), with 500 virtual patients per dose. Red diamonds are the published geometric mean anchors.

  Key takeaway: The boxes are nearly the same height across a 5.3-fold dose range. AUC hovers around 65–85 hr·µg/mL
  regardless of dose. Cmax is similarly compressed (~7.5–9 µg/mL). This is the foundational observation — dose
  escalation doesn't meaningfully increase exposure because of saturable absorption.

  ---
  Figure 2: saturable_bioavailability.png

  What it shows: Two panels explaining why exposure is flat.

  - Left panel: Relative bioavailability (F/F_180) vs dose. The curve drops steeply — at 960 mg the gut only absorbs
  ~20% as efficiently as at 180 mg. Each colored dot is a dose level.
  - Right panel: The consequence — AUCss vs dose. The solid blue line (saturable model) plateaus around 65–85
  hr·µg/mL. The dashed gray line shows what dose-proportional would look like (linear up to ~350 hr·µg/mL at 960 mg).
  Red diamonds = published anchors, blue dots = simulated geometric means.

  Key takeaway: The absorption machinery saturates at low doses. Giving 5x more drug barely changes what reaches the
  blood.

  ---
  Figure 3: er_efficacy_panel.png

  What it shows: Three panels on efficacy.

  - Panel A (Dose vs ORR): Paired bars — blue = published ORR, red = simulated ORR — at each dose. Error bars are
  binomial 95% CIs. N labels sit above each published bar. The 720 mg bar hits 50% but its CI spans nearly 0–80%
  (N=6). The 960 mg published (36%, N=124) has a tight CI.
  - Panel B (AUC Quartile ORR, 960 mg): ORR by AUC quartile within the 960 mg arm. Q1 (~22%) to Q4 (~34%) shows a
  shallow positive trend. Title says "model-implied, virtual patients."
  - Panel C (PFS by dose): Simple bar chart — published median PFS: 240 mg = 5.6 months, 960 mg = 5.4 months. Nearly
  identical.

  Key takeaway: ORR is ~30% across all doses (simulated). The published variability (25–50%) is noise from tiny sample
   sizes. PFS is essentially the same at 240 and 960 mg.

  ---
  Figure 4: dose_vs_exposure_vs_response.png

  What it shows: The three-link causal chain, one panel per link.

  - Panel A (Dose → Exposure): GM AUC vs dose with IQR shading. The line is nearly flat with a slight hump. Annotated
  "Saturable absorption." This is where the chain breaks.
  - Panel B (Exposure → Response): The Emax curve (red) showing ORR vs AUC. All five dose-level points cluster in a
  narrow AUC band (yellow shaded region, ~65–85 hr·µg/mL), sitting on the plateau of the curve. Even though a true E-R
   relationship exists, you can't see it because all doses give the same exposure.
  - Panel C (Dose → Response): The end result — a flat bar chart of simulated ORR by dose. The chain: flat input
  (Panel A) through a real curve (Panel B) produces a flat output (Panel C).

  Key takeaway: This is the teaching figure. The dose-response is flat not because the drug doesn't work, but because
  the pharmacokinetic step saturates. If you could somehow increase exposure, you would see more response — but dose
  escalation can't get you there.

  ---
  Figure 5: er_safety_panel.png

  What it shows: Three panels on hepatotoxicity.

  - Panel A (Cmax Quartile): Grade 3+ hepatotoxicity rate by Cmax quartile within the 960 mg arm. Clear gradient: Q1 =
   8%, Q2 = 12%, Q3 = 17%, Q4 = 27%. Color-coded green→red. Title says "model-implied, virtual patients."
  - Panel B (CPI Interaction): Four bars — 240 mg CPI−, 240 mg CPI+, 960 mg CPI−, 960 mg CPI+. Without prior CPI,
  rates are ~5–11%. With prior CPI, rates jump to ~21–30%. The dashed line is the published 16% pooled rate. Prior CPI
   is a much stronger risk factor than dose.
  - Panel C (Risk Stratification): A 2×2 matrix (low/high Cmax × CPI−/CPI+). Low Cmax + CPI− = 4.1%. High Cmax + CPI+
  = 41.7%. A 10-fold spread in risk based on two identifiable patient characteristics.

  Key takeaway: Hepatotoxicity has a real Cmax-driven exposure-safety relationship (unlike efficacy). Prior CPI
  dramatically amplifies risk. This creates an actionable risk stratification: patients with recent CPI and high Cmax
  are the ones most at risk.

  ---
  Figure 6: dose_optimization_dashboard.png

  What it shows: Six-panel dashboard comparing all dose levels.

  - A (ORR): ~30% across all doses. 180 mg is slightly higher (34%) by chance.
  - B (Hepatotox): 180 mg highest (~19%), others 15–16%. Higher Cmax at lower doses (saturable absorption) pushes the
  180 mg safety rate up.
  - C (Therapeutic Index): ORR / hepatotox rate. Ranges from 1.8 (180 mg) to 2.0 (960 mg). All similar.
  - D (Net Clinical Benefit): ORR minus hepatotox rate. All doses in the 13–15% range.
  - E (NNT vs NNH): NNT ~3.3 across all doses. NNH ~5.4–6.5. The gap between bars is the "margin of benefit."
  - F (Summary text): Head-to-head 240 vs 960: ORR 30% vs 30%, hepatotox 16% vs 15%, PFS 5.6 vs 5.4 mo. Conclusion:
  similar.

  Key takeaway: No dose level clearly dominates. The benefit-risk profile is approximately equivalent across the
  entire dose range, which is the core argument for why 960 mg is unnecessary.

  ---
  Figure 7: benefit_risk_overlay.png

  What it shows: A scatter plot with hepatotoxicity rate (x-axis) vs ORR (y-axis). Each dose is a colored circle with
  binomial 95% CI error bars on both axes. Published trial values are red/orange diamonds. Dashed diagonal lines mark
  constant "net clinical benefit" (NCB = ORR − hepatotox). The ideal corner is top-left (high efficacy, low toxicity).

  Key takeaway: All simulated dose points cluster tightly in the same region (~15–19% hepatotox, ~30–34% ORR). The
  published 240 mg diamond (24.8% ORR, 14% hepatotox) sits slightly below and left. The published 960 mg diamond
  (32.7% ORR, 16% hepatotox) is right among the simulated points. No dose clearly reaches the "ideal" corner. The
  clustering confirms there's no meaningful separation in benefit-risk across doses.

  ---
  Figure 8: covariate_forest_plot.png

  What it shows: A forest plot of AUC ratio vs reference (population GM) for 16 covariate subgroups at 960 mg. Each
  row is a horizontal line (90% CI) with a point estimate. Blue = CI includes 1.0. Red = CI excludes 1.0. The green
  band marks the 0.80–1.25 bioequivalence corridor. Numerical values are printed on the right.

  Significant (red): ECOG 1 (0.84), ECOG 2 (0.78), Low albumin (0.81), Body weight 100 kg (1.15), Hepatic mild
  impairment (1.10).

  Not significant (blue): Body weight 50/70/90 kg, ECOG 0, Normal albumin, Prior CPI, Renal function (all strata),
  Hepatic normal.

  Key takeaway: ECOG and albumin are the clinically meaningful covariates — sicker patients (higher ECOG, lower
  albumin) clear sotorasib faster, getting less exposure. Weight, renal function, and prior CPI don't drive PK. This
  matches the published popPK findings from Nagase et al. 2025.

  ---
  Tables

  virtual_patients.csv (2,500 rows)

  The core dataset. One row per virtual patient. Columns: patient ID, dose (mg), AUCss, Cmax,ss, Ctrough,ss, weight
  (kg), ECOG (0/1/2), prior CPI (0/1), albumin (g/dL). 500 patients at each of five dose levels. All downstream
  analyses read from this file.

  pk_summary_by_dose.csv

  Five rows (one per dose). Reports geometric mean AUC, geometric mean Cmax, median AUC, and CV% for AUC. Confirms the
   flat dose-exposure pattern: GM AUC ranges 65–86 hr·µg/mL, CV% around 70–81%.

  pk_anchor_vs_simulated.csv

  The PK calibration audit table. Five rows. For each dose: published GM AUC, simulated GM AUC, their ratio (all
  1.00–1.06), and the same for Cmax (all 0.99–1.03). This is the proof that the virtual population reproduces
  published PK.

  efficacy_anchor_vs_simulated.csv

  The efficacy calibration audit table. Five rows. Published ORR vs simulated ORR at each dose with ratio. The 180 mg
  and 960 mg ratios are 0.87–0.89 (close). The 240/360 mg ratios are ~1.3 (simulated higher than published) and 720 mg
   is 0.60 (simulated lower). This is expected — the published 720 mg ORR of 50% comes from just 3 out of 6 patients,
  and the Emax model can't reproduce that sampling noise.

  logistic_regression_orr.csv

  Two rows: intercept and AUC coefficient from the logistic regression (ORR ~ AUCss). AUC coefficient = 0.00182, p =
  0.0009, OR per 10 hr·µg/mL = 1.018. Statistically significant but clinically tiny — a 10-unit AUC increase predicts
  only 1.8% higher odds of response. This confirms the shallow E-R on the plateau.

  orr_by_exposure_quartile.csv

  Four rows (Q1–Q4) within the 960 mg arm. Median AUC ranges from 31 (Q1) to 130 (Q4). ORR: Q1 = 21.6%, Q2 = 36.8%, Q3
   = 32.0%, Q4 = 34.4%. Shallow positive trend — lowest quartile is lower, but Q2–Q4 are similar. All model-implied
  from the virtual population.

  pfs_descriptive.csv

  Two rows. Published median PFS: 960 mg = 5.4 months (N=104), 240 mg = 5.6 months (N=105). Implied HR = 1.037. No
  simulated data — these are just the published numbers passed through, confirming near-identical PFS.

  safety_anchor_vs_simulated.csv

  Five rows. Only 240 mg and 960 mg have published targets. 960 mg: published 16.0%, simulated 16.0%, ratio = 1.00
  (exact hit). 240 mg: published 14.0%, simulated 16.8%, ratio = 1.20 (at the acceptance boundary). Other doses show
  simulated rates only (no published reference).

  safety_by_exposure_quartile.csv

  Four rows (Q1–Q4) within the 960 mg arm by Cmax quartile. Median Cmax ranges from 4.1 (Q1) to 13.7 (Q4).
  Hepatotoxicity rate: Q1 = 8.0%, Q2 = 12.0%, Q3 = 16.8%, Q4 = 27.2%. Clear exposure-safety gradient — unlike
  efficacy, safety has a real Cmax-driven signal.

  safety_logistic_model.csv

  Three rows: intercept, Cmax coefficient, and prior CPI coefficient. Cmax OR = 1.08 per µg/mL (p < 0.0001). Prior CPI
   OR = 5.52 (p < 0.0001). Prior CPI is by far the dominant risk factor — it multiplies the odds of Grade 3+
  hepatotoxicity by 5.5x.

  dose_optimization_summary.csv

  Five rows (one per dose). The master benefit-risk table. Every dose shows ORR ~30%, hepatotox 13–19%, therapeutic
  index 1.8–2.0, NNT ~3.3, NNH ~5.4–6.5, net clinical benefit ~13–15%. No dose stands out. This is the quantitative
  backbone of the "240 mg suffices" argument.

  benefit_risk_comparison.csv

  Direct 240 mg vs 960 mg head-to-head. Simulated metrics plus published ORR and PFS for context. ORR: 29.8% vs 30.0%.
   Hepatotox: 16.2% vs 15.4%. Published PFS: 5.6 vs 5.4 months. The two doses are essentially interchangeable on every
   metric.

  covariate_effects.csv

  Sixteen rows (one per covariate subgroup). Columns: covariate name, subgroup, N, GM AUC, ratio vs reference, 90% CI
  bounds, significance flag. This is the data behind the forest plot. Five subgroups flagged "Yes" (CI excludes 1.0),
  two flagged "Ref" (reference categories), nine flagged "No."

