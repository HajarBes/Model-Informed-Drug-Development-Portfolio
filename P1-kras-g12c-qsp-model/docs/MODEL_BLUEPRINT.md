# MODEL BLUEPRINT — KRAS G12C QSP Model (M2-M6)

**Date:** 2026-02-12
**Architecture:** Single integrated ODE system (PK → KRAS → Signaling → Tumor)
**Implementation:** Python + scipy.integrate.solve_ivp

---

## Overview: State Variables (11 total)

| # | Variable | Description | Units | Module |
|---|----------|-------------|-------|--------|
| 1 | A_gut | Drug amount in gut | mg | PK (M1) |
| 2 | A_central | Drug amount in central compartment | mg | PK (M1) |
| 3 | KRAS_GDP | Inactive KRAS (drug-targetable) | µM | KRAS (M2) |
| 4 | KRAS_GTP | Active KRAS (signaling-competent) | µM | KRAS (M2) |
| 5 | KRAS_drug | Covalently bound KRAS-drug complex | µM | KRAS (M2) |
| 6 | R | Receptor drive (RTK/EGFR proxy) | dimensionless | Feedback (M3) |
| 7 | E | Pathway output (ERK proxy) | dimensionless | Feedback (M3) |
| 8 | T | Tumor volume | cm³ | Tumor (M4) |

Drug concentration derived: C(t) = A_central / V [mg/L], then converted to µM via MW.

---

## M2: KRAS GDP/GTP Cycling + Drug Binding

### Biology
KRAS cycles between inactive (GDP-bound, OFF) and active (GTP-bound, ON) states. GEFs promote GDP→GTP exchange; intrinsic GTPase + GAPs promote GTP→GDP hydrolysis. KRAS G12C has ~72% of WT GTPase activity (slower hydrolysis → more time in GTP state → oncogenic). Sotorasib binds covalently and irreversibly to KRAS-G12C in the GDP-bound state only.

### Equations

```
dKRAS_GDP/dt = k_hyd * KRAS_GTP                    # hydrolysis: GTP→GDP
             - k_GEF * R * KRAS_GDP                 # exchange: GDP→GTP (receptor-driven)
             - k_bind * C_drug_uM * KRAS_GDP         # covalent drug binding (irreversible)
             + k_syn                                  # protein synthesis
             - k_deg * KRAS_GDP                       # protein degradation

dKRAS_GTP/dt = k_GEF * R * KRAS_GTP_scaling * KRAS_GDP   # exchange: GDP→GTP
             - k_hyd * KRAS_GTP                            # hydrolysis: GTP→GDP
             - k_deg * KRAS_GTP                            # degradation

dKRAS_drug/dt = k_bind * C_drug_uM * KRAS_GDP       # covalent binding (irreversible)
              - k_deg * KRAS_drug                     # degraded along with protein turnover
```

### Parameters

| Parameter | Symbol | Value | Units | Source / Rationale |
|-----------|--------|-------|-------|--------------------|
| GTPase hydrolysis rate | k_hyd | 0.2 | /h | Stites 2018: intrinsic + GAP, ~72% of WT |
| GEF exchange rate (basal) | k_GEF | 0.1 | /h | Tuned so steady-state KRAS_GTP/total ~ 10-15% |
| Covalent binding rate | k_bind | 10.0 | /(µM·h) | From IC50 ~ 48 nM; tuned for ~90% TE at Cmax |
| Protein synthesis rate | k_syn | 0.029 | µM/h | k_deg × KRAS_total_ss = 0.029 × 1.0 |
| Protein degradation rate | k_deg | 0.029 | /h | t1/2 = 24h → ln(2)/24 |
| KRAS total (steady state) | KRAS_total_ss | 1.0 | µM | Stites 2018, Sumi 2021 |

### Target Engagement

```
TE(t) = KRAS_drug(t) / (KRAS_GDP(t) + KRAS_GTP(t) + KRAS_drug(t))
```

### Validation Checkpoint (M2)
- [ ] At steady state (no drug): KRAS_GTP ~ 10-15% of total
- [ ] After sotorasib 960 mg QD × 5 days: TE > 80% at Cmax, oscillates with dosing
- [ ] KRAS_GTP suppressed to < 30% of baseline at peak drug exposure

### Expected Figure: Fig 2 — Target Engagement
- Top panel: Drug concentration C(t) over 5 days
- Bottom panel: TE(t) and KRAS_GTP(t)/KRAS_GTP_baseline showing oscillating suppression

---

## M3: Receptor Feedback + ERK Rebound

### Biology
ERK (pathway output) normally suppresses RTK signaling via negative feedback (DUSP, SPRY, etc.). When KRAS is inhibited → ERK drops → feedback released → RTK drive increases → wild-type RAS activates → ERK rebounds. This is the core adaptive resistance mechanism.

We model this as a 2-variable feedback loop:
- **R(t)**: receptor drive (RTK/EGFR proxy). Increases when E is low (feedback release).
- **E(t)**: pathway output (ERK proxy). Driven by KRAS_GTP and R.

### Equations

```
dR/dt = (1/tau_R) * (R_basal + G_fb * (1 - E/E_ss) - R)
```
Interpretation: R relaxes toward `R_basal + G_fb * (1 - E/E_ss)` with time constant tau_R.
When E drops below steady state → feedback term positive → R increases.

```
dE/dt = (1/tau_E) * (R * (KRAS_GTP / KRAS_GTP_ss) - E)
```
Interpretation: E is driven by R × (normalized KRAS_GTP). Relaxes with time constant tau_E.

At steady state (no drug): R = R_basal (E = E_ss, feedback term = 0), and E = R_basal × 1.0 = E_ss.
Choose R_basal = E_ss = 1.0 (dimensionless, normalized).

### Parameters

| Parameter | Symbol | NSCLC | CRC | Units | Rationale |
|-----------|--------|-------|-----|-------|-----------|
| Receptor time constant | tau_R | 12 | 12 | h | RTK upregulation timescale |
| Pathway output time constant | tau_E | 2 | 2 | h | ERK signaling is fast |
| Basal receptor drive | R_basal | 1.0 | 1.0 | — | Normalized |
| **Feedback gain** | **G_fb** | **0.5** | **2.0** | — | **Key lineage parameter**: CRC has 4× stronger feedback than NSCLC |
| Steady-state pathway output | E_ss | 1.0 | 1.0 | — | Normalized |
| Steady-state KRAS_GTP | KRAS_GTP_ss | set from M2 | set from M2 | µM | From M2 no-drug steady state |

### EGFR Blocker Effect (M6 preview)
When panitumumab is present, receptor drive is attenuated:
```
R_effective = R * (1 - f_block)
```
where f_block ∈ [0, 1] represents EGFR blockade efficacy.
For the portfolio model, we use a simple step: f_block = 0.6 when panitumumab is dosed.

### Validation Checkpoint (M3)
- [ ] No drug: R = 1.0, E = 1.0 at steady state
- [ ] After KRAS inhibition: E drops then rebounds toward ~75% baseline by 72h (matches Xue 2023)
- [ ] CRC rebound is stronger than NSCLC (G_fb = 2.0 vs 0.5)
- [ ] R increases when E drops (feedback release visible in R(t) plot)

### Expected Figure: Fig 3 — Adaptive Rebound
- Panel A: E(t) for NSCLC mono — moderate rebound
- Panel B: E(t) for CRC mono — strong rebound to ~75-80% baseline
- Panel C: R(t) showing receptor drive surge in CRC > NSCLC

---

## M4: Tumor Growth Module

### Biology
Tumor growth depends on pathway output E(t). Higher E → more proliferation. The model uses a logistic growth framework where the net growth rate is modulated by E(t).

### Equation

```
dT/dt = rho_grow * (E / E_ss)^n_hill * T * (1 - T/T_max) - rho_kill * max(0, 1 - E/E_ss) * T
```

- **Growth term:** logistic growth with rate modulated by E. When E = E_ss (baseline), growth rate = rho_grow.
- **Kill term:** activated only when E < E_ss (drug effect). Stronger suppression → more kill.
- **n_hill:** Hill coefficient for E-driven proliferation (nonlinearity)

### Parameters

| Parameter | Symbol | Value | Units | Rationale |
|-----------|--------|-------|-------|-----------|
| Baseline growth rate | rho_grow | 0.005 | /h (~doubling time 6 days) | Typical solid tumor |
| Max kill rate | rho_kill | 0.003 | /h | Tuned for ORR-consistent shrinkage |
| Hill coefficient | n_hill | 2 | — | Moderate nonlinearity |
| Carrying capacity | T_max | 100 | cm³ | Arbitrary normalization |
| Initial tumor volume | T_0 | 10 | cm³ | Starting tumor (assessable disease) |
| Pathway output baseline | E_ss | 1.0 | — | Normalized |

### Validation Checkpoint (M4)
- [ ] No drug: tumor grows logistically toward T_max
- [ ] NSCLC mono: tumor shrinks initially, then regrows (partial response, then progression)
- [ ] Time to nadir ~2-3 months, consistent with PFS ~6.8 mo
- [ ] Best tumor change ~ -30 to -40% (consistent with 37% ORR)

### Expected Figure: Fig 4 — Tumor Growth Curves
- Vehicle (no drug): monotonic growth
- NSCLC mono: initial shrinkage, nadir, regrowth
- CRC mono: minimal shrinkage, faster regrowth (weaker response)

---

## M5: CRC vs NSCLC Parameterization

### Key Difference: Feedback Gain

| Parameter | NSCLC | CRC | Biological basis |
|-----------|-------|-----|-----------------|
| G_fb | 0.5 | 2.0 | CRC has stronger RTK/EGFR-driven feedback |
| rho_kill | 0.003 | 0.003 | Same drug effect on pathway |
| rho_grow | 0.005 | 0.006 | CRC slightly more aggressive (optional) |

All other parameters identical. The lineage difference is encoded entirely in G_fb (and optionally rho_grow), which is biologically justified:
- CRC has higher baseline EGFR expression
- Multiple RTKs contribute to feedback in CRC
- Clinical ORR: 37% (NSCLC) vs 10% (CRC) — ~4× difference

### Validation Checkpoint (M5)
- [ ] NSCLC mono: meaningful shrinkage in many virtual patients
- [ ] CRC mono: minimal/no shrinkage in most patients
- [ ] Qualitative ORR ratio roughly consistent with 37%:10%

### Expected Figure: Fig 5 — CRC vs NSCLC
- Side-by-side tumor curves: NSCLC mono vs CRC mono
- Annotation showing differential rebound (E(t) overlay)

---

## M6: Combination Therapy (EGFR Blockade)

### Implementation

Add f_block parameter to the receptor drive equation:

```
dR/dt = (1/tau_R) * (R_basal + G_fb * (1 - E/E_ss) - R) * (1 - f_block)
```

Or equivalently, the R set-point is reduced:
```
R_target = (R_basal + G_fb * (1 - E/E_ss)) * (1 - f_block)
dR/dt = (1/tau_R) * (R_target - R)
```

**Panitumumab dosing:** 6 mg/kg IV q2w. For simplicity, model as constant f_block when drug is present (panitumumab t1/2 ~7.5 days → relatively flat trough with q2w dosing).

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| f_block | 0.6 | EGFR blockade reduces receptor drive by 60%. Tuned for ORR ~26% in CRC combo. |

### Validation Checkpoint (M6)
- [ ] CRC combo: better tumor shrinkage than CRC mono
- [ ] ERK rebound suppressed in combo vs mono
- [ ] ORR qualitatively consistent: ~10% (mono) → ~26% (combo)

### Expected Figure: Fig 6 — Combination Benefit
- CRC mono vs CRC combo tumor curves
- Inset: E(t) showing rebound suppression in combo

---

## M7-M8: Virtual Population + Sensitivity (Optional)

### M7: Virtual Population
- Sample G_fb, rho_grow, rho_kill, k_bind from log-normal distributions
- N = 500-2000 virtual patients
- Define "responder" as max tumor shrinkage > 30% (RECIST partial response proxy)
- Generate waterfall plot of best tumor change

### M8: Global Sensitivity (Morris or Sobol)
- Outcome: minimum tumor volume, time-to-regrowth, ERK rebound amplitude
- Parameters: k_GEF, G_fb, k_bind, rho_grow, rho_kill, tau_R
- Generate tornado/ranking plot

---

## Implementation Architecture

### Single file: `src/kras_qsp.py`

```python
# Pseudocode structure

@dataclass
class Params:
    # PK
    dose_mg, tau_h, ka, ke, V, F, MW
    # KRAS
    k_hyd, k_GEF, k_bind, k_syn, k_deg, KRAS_total_ss
    # Feedback
    tau_R, tau_E, R_basal, G_fb, E_ss
    # Tumor
    rho_grow, rho_kill, n_hill, T_max, T_0
    # Combo
    f_block

def rhs(t, y, params):
    """ODE right-hand side for full integrated system."""
    A_gut, A_cent, KRAS_GDP, KRAS_GTP, KRAS_drug, R, E, T = y

    # PK
    C_mgL = A_cent / params.V
    C_uM = C_mgL * 1000 / params.MW  # mg/L → µM

    # Dosing events handled externally (events or piecewise)

    # KRAS cycling
    dKRAS_GDP = (k_hyd * KRAS_GTP
                 - k_GEF * R * KRAS_GDP
                 - k_bind * C_uM * KRAS_GDP
                 + k_syn - k_deg * KRAS_GDP)
    dKRAS_GTP = (k_GEF * R * KRAS_GDP
                 - k_hyd * KRAS_GTP
                 - k_deg * KRAS_GTP)
    dKRAS_drug = (k_bind * C_uM * KRAS_GDP
                  - k_deg * KRAS_drug)

    # Feedback
    KRAS_GTP_ss = ...  # from initial steady state
    R_target = (R_basal + G_fb * max(0, 1 - E/E_ss)) * (1 - f_block)
    dR = (1/tau_R) * (R_target - R)
    dE = (1/tau_E) * (R * (KRAS_GTP / KRAS_GTP_ss) - E)

    # Tumor
    growth = rho_grow * (E/E_ss)**n_hill * T * (1 - T/T_max)
    kill = rho_kill * max(0, 1 - E/E_ss) * T
    dT = growth - kill

    return [dA_gut, dA_cent, dKRAS_GDP, dKRAS_GTP, dKRAS_drug, dR, dE, dT]

def simulate(params, scenario, t_days=180):
    """Run simulation with dosing events."""
    # ...solve_ivp with event handling for doses...

def plot_all(results, scenario_name, outdir):
    """Generate publication-quality figures."""
    # ...7 figure functions...
```

### CLI Interface

```bash
python src/kras_qsp.py --scenario nsclc_mono    # Fig 1-4
python src/kras_qsp.py --scenario crc_mono      # Fig 4-5
python src/kras_qsp.py --scenario crc_combo     # Fig 6
python src/kras_qsp.py --scenario all           # All scenarios + comparison
python src/kras_qsp.py --vpop --n 500           # Fig 7 (M7)
python src/kras_qsp.py --sensitivity            # Fig 8 (M8)
```

---

## Figure Plan (7 Publication-Quality Figures)

| Fig | Title | Content | Module |
|-----|-------|---------|--------|
| 1 | PK Profile | C(t) over 7 days, steady-state profile | M1 |
| 2 | Target Engagement | TE(t), KRAS_GTP(t) over 7 days | M2 |
| 3 | Adaptive Rebound | E(t), R(t) for NSCLC vs CRC | M3 |
| 4 | Tumor Growth Curves | Vehicle, NSCLC mono, CRC mono | M4 |
| 5 | CRC vs NSCLC Comparison | Side-by-side: tumor + E(t) + R(t) | M5 |
| 6 | Combination Benefit | CRC mono vs combo: tumor + rebound | M6 |
| 7 | Model Schematic | Pathway diagram with state variables | — |

---

## Build Order

1. Implement integrated ODE (PK + KRAS + feedback + tumor) in single `kras_qsp.py`
2. Verify no-drug steady state (M2 + M3 checkpoints)
3. Run NSCLC mono scenario → validate TE, rebound, tumor (Fig 1-4)
4. Run CRC mono → validate differential response (Fig 5)
5. Run CRC combo → validate rebound suppression and tumor improvement (Fig 6)
6. Generate all figures
7. Update `model/equations.md` with final equations
