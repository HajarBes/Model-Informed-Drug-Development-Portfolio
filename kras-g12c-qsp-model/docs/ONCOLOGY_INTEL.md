# ONCOLOGY INTEL — KRAS G12C QSP Literature & Regulatory Review

**Date:** 2026-02-12
**Purpose:** Inform model structure, parameterization, and regulatory alignment for KRAS-G12C QSP portfolio model.

---

## 1. Key Papers Summary

### 1a. QSP / Mechanistic Modeling Papers

| # | Paper | Year | Journal | Relevance | Key Contribution |
|---|-------|------|---------|-----------|-----------------|
| 1 | **Stites & Shaw** — "Quantitative Systems Pharmacology Analysis of KRAS G12C Covalent Inhibitors" | 2018 | CPT: Pharmacometrics & Systems Pharmacology | **Foundational** | First QSP model of KRAS GDP/GTP cycling with covalent inhibitor binding. Includes protein turnover, NPI vs SIIPI classification. Parameters in Table 1. |
| 2 | **Sumi et al.** — "Virtual clinical trial simulations for a novel KRASG12C inhibitor (ASP2453) in NSCLC" | 2021 | PMC8376128 | **Directly relevant** | Full QSP: RAS cycling → signal transduction (RAF/MEK/ERK/S6 + PI3K/AKT) → tumor growth. 27 species, 75 parameters. Virtual population of 1000 patients. Clinical validation against CodeBreaK 100 waterfall plots. |
| 3 | **Xue et al.** — "KRASG12C-independent feedback activation of wild-type RAS constrains KRASG12C inhibitor efficacy" | 2023 | Cell Reports | **Critical for M3** | Wild-type NRAS/HRAS reactivation (not mutant KRAS shift) drives ERK rebound to ~75% of baseline by 72h. Multiple RTKs mediate feedback — not EGFR alone. |

### 1b. Biology / Resistance Mechanism Papers

| # | Paper | Year | Journal | Key Findings |
|---|-------|------|---------|-------------|
| 4 | **Amodio et al.** — "Vertical pathway inhibition overcomes adaptive feedback resistance to KRASG12C inhibition" | 2020 | Clin Cancer Res | ERK rebound via RTK→wild-type RAS. SHP2 inhibitor blocks RTK-mediated WT-RAS activation. Vertical inhibition (SHP2i or MEKi) more durable than horizontal. |
| 5 | **Ryan & Corcoran** — "EGFR-inhibition in CRC KRAS G12C models" | 2021 | Clin Cancer Res | Panitumumab effective only in CRC subset (EGFR-dependent feedback). ~1/3 CRC models don't respond to EGFR blockade alone. Cell-type-specific adaptive signaling. |
| 6 | **Riedl et al.** — "Emerging landscape of KRAS inhibitors in cancer treatment" | 2026 | Cancer Cell (review) | Comprehensive review: KRAS biology, inhibitor classes (OFF-state vs ON-state), resistance mechanisms, clinical outcomes across tumor types, combination strategies. **Already in repo as `mmc2 (2).pdf`.** |

### 1c. Clinical Data Papers

| # | Paper | Year | Trial | Indication | Key Data |
|---|-------|------|-------|-----------|---------|
| 7 | **Hong et al.** — CodeBreaK 100 (NSCLC) | 2021 | NEJM | NSCLC mono | **ORR 37.1%**, mPFS 6.8 mo |
| 8 | **Fakih et al.** — CodeBreaK 100 (CRC cohort) | 2022 | Lancet Oncol | **CRC mono** | **ORR 9.7%** (6/62 pts), mPFS 4.0 mo |
| 9 | **Fakih et al.** — CodeBreaK 300 | 2023 | NEJM | CRC combo | **ORR 26.4%**, mPFS 5.6 mo (soto + pani vs standard) |
| 10 | **Canon et al.** — "Discovery of AMG 510" | 2019 | Nature | Preclinical | First disclosure of sotorasib. Xenograft data. IC50 ~48 nM (H358 cells). |

---

## 2. Extracted Model Parameters

### 2a. KRAS Biochemistry (from Stites 2018 + literature)

| Parameter | Symbol | Value | Units | Source |
|-----------|--------|-------|-------|--------|
| KRAS total concentration | [KRAS]_total | ~1 × 10⁶ | molecules/cell (~1 µM) | Stites 2018, Sumi 2021 |
| GTP hydrolysis rate (WT KRAS) | k_GAP | ~0.02 | /s | Stites 2018 |
| GTP hydrolysis rate (G12C) | k_GAP,G12C | 72% of WT | relative | Stites 2018 (ScalingfactorGTPase) |
| GEF exchange rate | k_GEF | ~0.001-0.01 | /s | Literature range |
| Cellular GTP concentration | [GTP] | 180 µM | M | Stites 2018 |
| Cellular GDP concentration | [GDP] | 18 µM | M | Stites 2018 |
| KRAS protein half-life | t1/2_KRAS | 24 | h | Sumi 2021 |
| Ras degradation rate | k_deg | 8 × 10⁻⁶ | /s (~0.029/h) | Stites 2018 |

### 2b. Drug Parameters (Sotorasib)

| Parameter | Symbol | Value | Units | Source |
|-----------|--------|-------|-------|--------|
| PK half-life | t1/2 | 5 | h | FDA label |
| Dose | — | 960 mg | QD oral | FDA label |
| Covalent binding rate (SIIPI class) | k_on | 76 | /Ms (ARS-853 proxy) | Stites 2018 |
| Cellular IC50 (pERK, H358) | IC50 | 48 nM | nM | Canon 2019 |
| MW | — | 560.6 | g/mol | PubChem |
| Vss/F | — | 211 | L | FDA label |
| Bioavailability | F | ~0.73 | — | FDA label |

### 2c. Signaling / Feedback Parameters (from Sumi 2021)

| Parameter | Description | Value | Source |
|-----------|-------------|-------|--------|
| ERK rebound to baseline | pERK recovery after 72h treatment | ~75% of baseline | Xue 2023 |
| WT-RAS activation fold-change | HRAS/NRAS-GTP increase upon KRAS-G12C inhibition | ~3-fold | Xue 2023 |
| Feedback loops in model | DUSP, SPRY, cMYC, AKT→cRAF | 4 loops | Sumi 2021 |

### 2d. Tumor Parameters

| Parameter | Description | Value | Source |
|-----------|-------------|-------|--------|
| NSCLC ORR (mono) | Sotorasib 960 mg | 37.1% | CodeBreaK 100 |
| CRC ORR (mono) | Sotorasib 960 mg | 9.7% | CodeBreaK 100 CRC |
| CRC ORR (combo) | Sotorasib + panitumumab | 26.4% | CodeBreaK 300 |
| CRC mPFS (combo) | Sotorasib + panitumumab | 5.6 mo | CodeBreaK 300 |
| NSCLC mPFS (mono) | Sotorasib 960 mg | 6.8 mo | CodeBreaK 100 |

---

## 3. Model Structure Benchmarking

### Published QSP architectures for KRAS G12C:

**Stites 2018 (foundational):**
```
[Drug] --kon--> KRAS_GDP:Drug (covalent, irreversible)
KRAS_GDP <--kGEF/kGAP--> KRAS_GTP --> Effector binding
+ protein turnover (synthesis/degradation)
```
- 16 species, mass-action kinetics
- No downstream signaling or tumor module
- Focus on target engagement and inhibitor class comparison

**Sumi 2021 (full QSP):**
```
PK (oral) --> [Drug] --> KRAS cycling --> RAF/MEK/ERK/S6 --> Tumor growth
                                      --> PI3K/AKT -----------^
              + 4 negative feedback loops (DUSP, SPRY, cMYC, AKT→cRAF)
```
- 27 species, 75 parameters
- Logic-based signal transduction (not full ODE)
- Tumor growth: proliferation driven by pS6, death by carrying capacity
- Virtual population: 1000 patients, prevalence-weighted

### Our model positioning:

| Feature | Stites 2018 | Sumi 2021 | **Our Model (planned)** |
|---------|-------------|-----------|------------------------|
| PK | No | Yes (1-comp) | Yes (1-comp oral) |
| KRAS cycling | Full mass-action | Simplified | ODE mass-action |
| Drug binding | Covalent (irreversible) | Covalent | Covalent to GDP-bound |
| Downstream signaling | No | Logic-based (4 feedbacks) | ODE-based (simplified: R→E feedback) |
| Tumor growth | No | Yes (pS6-driven) | Yes (E-driven) |
| Lineage comparison | No | No (NSCLC only) | **Yes (CRC vs NSCLC)** |
| Combo therapy | No | No | **Yes (EGFR blockade)** |
| Virtual population | No | Yes (1000 pts) | Optional (M7) |
| Regulatory alignment | Research | Drug development | **Portfolio / Project Optimus** |

**Our competitive advantage:** No published QSP model simultaneously captures CRC vs NSCLC lineage differences AND combination benefit with EGFR blockade. This is the niche.

---

## 4. Regulatory Checklist for KRAS Combo Submissions

### 4a. FDA Project Optimus Alignment

Project Optimus (final guidance Aug 2024) requires:
- [ ] Dose-exposure-response characterization (not just MTD)
- [ ] Randomized dose comparison in Phase II
- [ ] Model-informed dose selection justification
- [ ] QSP/PBPK models accepted as supporting evidence for dose selection

**Relevance to our model:** A QSP model linking PK exposure → target engagement → pathway rebound → tumor response directly supports Project Optimus dose-response characterization. Showing that combo (soto + pani) suppresses rebound and improves tumor control at the same sotorasib dose is a compelling regulatory narrative.

### 4b. KRAS Inhibitor Regulatory Milestones

| Date | Event | Relevance |
|------|-------|-----------|
| May 2021 | Sotorasib accelerated approval (NSCLC) | First KRAS-targeted therapy |
| Dec 2022 | Adagrasib accelerated approval (NSCLC) | Second KRAS G12C agent |
| Jan 2025 | **Sotorasib + panitumumab approved (mCRC)** | First KRAS combo approval; CodeBreaK 300 |
| 2025-2026 | Next-gen inhibitors (divarasib, RAS-ON inhibitors) | Expanding landscape |

### 4c. QSP Model Credibility Expectations (FDA/EMA)

Neither FDA nor EMA has published KRAS-specific QSP guidance, but general expectations from MIDD framework:

1. **Context of Use (CoU):** Clearly state what the model is and isn't predicting
2. **Model qualification:**
   - Verification: code reproduces equations correctly
   - Validation: model outputs match independent data (e.g., ORR, PFS)
   - Sensitivity analysis: identify key drivers of uncertainty
3. **Parameter justification:** Every parameter traceable to published source
4. **Virtual population:** If used, describe sampling method and plausibility filtering
5. **Transparency:** Code available, equations documented, assumptions stated

### 4d. Credibility Checklist for Our Model

- [ ] All parameters sourced from published literature with citations
- [ ] PK validated against FDA label Cmax/AUC
- [ ] Target engagement consistent with IC50 data
- [ ] ERK rebound matches published ~75% recovery at 72h
- [ ] CRC mono ORR ~10% vs NSCLC mono ORR ~37% reproduced qualitatively
- [ ] Combo (CRC + EGFR blockade) shows improved response vs mono
- [ ] Sensitivity analysis identifies key drivers (k_GEF, feedback gain)
- [ ] Code documented and reproducible

---

## 5. Key Biological Insights for Model Design

### 5a. Why CRC responds worse than NSCLC to KRAS G12C monotherapy

1. **Stronger receptor-driven feedback in CRC:** CRC cells have higher baseline EGFR signaling. When KRAS G12C is inhibited, ERK drops → negative feedback on RTKs is released → RTK (especially EGFR) signaling surges → wild-type RAS (NRAS, HRAS) activates → ERK rebounds.
2. **Multiple RTKs involved:** Not just EGFR — HER2, HER3, FGFR, MET can all contribute. But EGFR is dominant in CRC.
3. **WT-RAS bypasses the drug:** Sotorasib only binds KRAS-G12C (GDP-bound). When WT-RAS is activated by RTK feedback, the drug has no target.

### 5b. Why combo (sotorasib + panitumumab) works in CRC

1. Panitumumab blocks EGFR → reduces the RTK-driven feedback
2. Less WT-RAS activation → less ERK rebound
3. More sustained pathway suppression → better tumor control
4. ORR jumps from ~10% (mono) to ~26% (combo)

### 5c. Key modeling implication

The model MUST encode:
- **Receptor drive R(t)** that increases when pathway output E(t) decreases (negative feedback)
- **Lineage parameter** controlling feedback gain (high in CRC, low in NSCLC)
- **EGFR blocker effect** that directly reduces R(t) in CRC

---

## Sources

- [Stites & Shaw 2018 — QSP Analysis of KRAS G12C Covalent Inhibitors](https://pmc.ncbi.nlm.nih.gov/articles/PMC5980551/)
- [Sumi et al. 2021 — Virtual Clinical Trial for ASP2453 in NSCLC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8376128/)
- [Xue et al. 2023 — WT-RAS Feedback Constrains KRASG12C Inhibitor Efficacy](https://pmc.ncbi.nlm.nih.gov/articles/PMC9809542/)
- [Amodio et al. 2020 — Vertical Pathway Inhibition Overcomes Adaptive Resistance](https://pmc.ncbi.nlm.nih.gov/articles/PMC7124991/)
- [Riedl et al. 2026 — Emerging Landscape of KRAS Inhibitors (Cancer Cell review)](https://doi.org/10.1016/j.ccell.2026.01.001)
- [Canon et al. 2019 — Discovery of AMG 510 (Sotorasib)](https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b01180)
- [Hong et al. 2021 — CodeBreaK 100 NSCLC (NEJM)](https://www.nejm.org/doi/full/10.1056/NEJMoa1917239)
- [Fakih et al. 2022 — CodeBreaK 100 CRC cohort (Lancet Oncol)](https://pubmed.ncbi.nlm.nih.gov/34919824/)
- [Fakih et al. 2023 — CodeBreaK 300 CRC Combo (NEJM)](https://www.nejm.org/doi/full/10.1056/NEJMoa2308795)
- [FDA Project Optimus — Oncology Dose Optimization](https://www.fda.gov/media/164555/download)
- [Gao et al. 2024 — Realizing Project Optimus (CPT:PSP)](https://ascpt.onlinelibrary.wiley.com/doi/10.1002/psp4.13079)
- [FDA — QSP/MIDD Industry Perspective 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11576823/)
