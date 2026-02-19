# Fa Definition — Fraction Absorbed from GI Lumen

## Definition

**Fa** is the fraction of the administered dose absorbed through the intestinal wall into the portal circulation (before first-pass metabolism).

$$
Fa(t) = \frac{1}{Dose} \int_0^t \sum_{i \in \text{absorbing}} k_{abs,i} \cdot D_i(\tau) \, d\tau
$$

where:
- $k_{abs,i} = P_{eff} \times SA_i / V_{lumen,i}$ (absorption rate constant for segment $i$, h$^{-1}$)
- $D_i(\tau)$ = dissolved drug mass in segment $i$ at time $\tau$ (mg)
- The sum is over absorbing segments only (duodenum through colon; stomach excluded)

## Computation

Fa is tracked as a **cumulative ODE state** (state index 16, the 17th state):

```
dy[16]/dt = sum_i( k_abs_i * D_i )
```

This integral accumulates mass absorbed through the intestinal wall over time. It is then normalized:

```
Fa = y[16] / dose
```

## Why not `1 - (S + D) / dose`?

The naive formula `1 - (S_total + D_total) / dose` computes the fraction of drug that **left the GI tract by any route**, including:
- Absorption through the intestinal wall (correct)
- **Fecal transit** out of the colon (incorrect — this is NOT absorption)

At late times (t >> transit time), both S and D approach zero as drug either absorbs or transits fecally. The naive formula converges to 1.0, overstating Fa.

## Consistency check

The regional absorption breakdown (`absorption_by_segment_fasted.csv`) computes absorbed mass per segment by integrating `k_abs_i * D_i * dt` post-hoc. The sum across all segments should equal `Fa * dose` from the ODE state.

| Quantity | Source | Value (10 mg, fasted, 48 h) |
|----------|--------|----------------------------|
| Fa (ODE state 17) | `cum_absorbed_mg / dose` | 0.9715 |
| Regional sum | `absorption_by_segment_fasted.csv` | 9.71 mg = 97.1% |
| Fecal loss | `dose - GI_remaining - absorbed` | 0.285 mg (2.85%) |
| F_oral | Fa * Fg * Fh = 0.971 * 0.50 * 0.50 | 0.243 |

## Related quantities

| Quantity | Definition | Formula |
|----------|------------|---------|
| **Fa** | Fraction absorbed through intestinal wall | ODE state 16 / dose |
| **Fg** | Gut wall availability (fraction surviving gut CYP3A4) | 0.50 (fixed, Lundahl 1997) |
| **Fh** | Hepatic availability (fraction surviving liver first-pass) | 0.50 (fixed, Lundahl 1997) |
| **F_oral** | Overall oral bioavailability | Fa * Fg * Fh |
