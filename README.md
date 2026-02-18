# Model-Informed Drug Development Portfolio

**Hajar Besbassi** | Pharmacometrics & QSP

This repository showcases quantitative systems pharmacology (QSP) and model-informed drug development (MIDD) projects built to industry and regulatory standards.

---

## Projects

### [KRAS G12C QSP Model](01-KRAS-G12C-qsp-model/)

A 9-state ODE model explaining why sotorasib (LUMAKRAS) achieves durable responses in NSCLC but not in CRC, and how EGFR blockade (panitumumab) rescues CRC outcomes through combination therapy.

**Highlights:**
- Mechanistic resistance framework with adaptive feedback rebound and bypass signaling
- Virtual population (N=200) with median PFS matching clinical data within 12%
- Morris global sensitivity analysis identifying key resistance drivers
- ICH M15-aligned credibility assessment and regulatory-style documentation

**Tools:** Python, SciPy, SALib | **Domain:** Oncology, Targeted Therapy, Adaptive Resistance

---

### [POPPK — Oral Midazolam CYP3A4 Probe](02-POPPK-midazolam-cyp3a4-probe/)

Population pharmacokinetic analysis of the standard CYP3A4 probe substrate (7.5 mg, N=120), with parameters calibrated to the OSP Midazolam PBPK model (qualified against multiple published clinical studies).

**Highlights:**
- 2-compartment SAEM estimation (nlmixr2/rxode2), dAIC = 467 vs 1-compartment
- Allometric weight scaling with verified exposure analysis and bootstrap CIs
- Production diagnostics: GOF, VPC, eta distributions, individual fits, forest plot
- CYP3A4 inhibition bridge connecting baseline PK to DDI risk
- External literature qualification against 12 OSP studies (11/12 within 2-fold AUC ratio)

**Tools:** R, nlmixr2, rxode2, Python | **Domain:** Clinical Pharmacology, PopPK, DDI

---

### [Static DDI Risk Assessment](03-static-DDI-framework/)

Drug-drug interaction screening per FDA 2020 guidance and ICH M12 (2024), with two cases: ketoconazole benchmark validation and sotorasib multi-mechanism perpetrator analysis.

**Highlights:**
- Reversible inhibition, TDI, and induction pathways with net effect logic
- Monte Carlo uncertainty quantification and tornado sensitivity analysis
- Identifies when PBPK is essential to resolve opposing mechanisms
- Regulatory-style case reports with full mathematical derivations

**Tools:** R | **Domain:** Clinical Pharmacology, Drug-Drug Interactions, Regulatory Science
