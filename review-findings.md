## REVIEW STATUS: 1 real P0 fixed, 1 real P1 fixed, 3 false positives identified
## Multi-Persona Review: adaptsim.html
### Date: 2026-03-24
### Summary: 3 P0 (1 real + 2 false positive), 5 P1 (2 real + 3 false positive), 4 P2

---

#### P0 -- Critical

- **P0-1** [Statistical]: OBF spending function uses **two-sided** alpha spending even for **one-sided** tests (line 1043-1044)
  - `alphaSpendOBF` computes `normalQuantile(1 - alpha/2)` regardless of sidedness. For one-sided alpha=0.025, it should use `normalQuantile(1 - alpha)` = `normalQuantile(0.975)`, not `normalQuantile(1 - 0.0125)` = `normalQuantile(0.9875)`.
  - The Pocock and Lan-DeMets functions are scale-free (just `alpha * f(t)`) so they're unaffected, but OBF's formula explicitly uses the quantile.
  - Impact: OBF boundaries are too conservative for one-sided tests (z1 too high at first look).
  - Fix: Pass `twoSided` flag to `alphaSpendOBF` and use `normalQuantile(1 - alpha)` for one-sided.

- **P0-2** [Statistical]: Monte Carlo drift parameterization inconsistency (line 1842 vs 1848-1849)
  - At look 1 (line 1842): `Z[0] = theta * sqrt(t_0) + normalRandom * sqrt(t_0)`
  - At look k>0 (line 1848-1849): `drift = theta * (t_k - t_prev) / sqrt(t_k)`, then `Z[k] = Z[k-1] * rho + drift + normalRandom * sigma`
  - The drift formula at k>0 is: `theta * dt / sqrt(t_k)`. Let's verify: E[Z_k] should be `theta * sqrt(t_k)`.
  - E[Z_k] = E[Z_{k-1}] * rho + theta * dt / sqrt(t_k) = theta*sqrt(t_{k-1}) * sqrt(t_{k-1}/t_k) + theta*dt/sqrt(t_k) = theta*t_{k-1}/sqrt(t_k) + theta*dt/sqrt(t_k) = theta*(t_{k-1}+dt)/sqrt(t_k) = theta*t_k/sqrt(t_k) = theta*sqrt(t_k). CORRECT.
  - At look 1: E[Z_0] = theta*sqrt(t_0). Var[Z_0] = t_0. For the Z-statistic to have unit variance under H0, we need Var[Z_0] = 1, not t_0. But `normalRandom * sqrt(t_0)` gives variance t_0.
  - **This means the variance of Z at look 1 is t_0, not 1.** For equal spacing with K=2, t_0=0.5, so Var[Z_0]=0.5 instead of 1. This makes early-look Z-stats have smaller variance → inflated Type I error at interim.
  - Fix: Z[0] should be `theta * sqrt(t_0) + normalRandom(rng)` (unit variance), then subsequent Z[k] = Z[k-1]*rho + drift + normalRandom*sigma.
  - Actually wait — let me re-check. The standard formulation has Z_k = B(t_k)/sqrt(t_k) where B is Brownian motion. Then Var(Z_k) = t_k/t_k = 1. The increment W_k = B(t_k) - B(t_{k-1}) ~ N(0, dt). Z_k = (sum W_j)/sqrt(t_k). So Z_1 = W_1/sqrt(t_1), Var(Z_1) = dt_1/t_1 = 1 (since dt_1 = t_1). The code has Z[0] = theta*sqrt(t_0) + normalRandom*sqrt(t_0). Since normalRandom ~ N(0,1), this gives Var = t_0. If t_0 = t_1 (first info frac), this equals t_1 which should be...
  - Actually: Var(Z_1) = Var(W_1/sqrt(t_1)) = t_1/t_1 = 1. The code gives Var = t_0 = t_1. So if t_1 = 0.5, Var = 0.5 ≠ 1. **This IS a bug.**
  - Suggested fix: `Z[0] = theta * Math.sqrt(t_k) + normalRandom(rng);` (remove the `* Math.sqrt(t_k)` on the noise term)

- **P0-3** [Security]: Web Worker Blob URL not revoked (lines 1955-1984)
  - `runSimulation` creates a Blob URL via `URL.createObjectURL(blob)` but never calls `URL.revokeObjectURL()`. Each simulation run leaks a Blob URL. After 20+ runs, this becomes a memory leak.
  - Fix: Store the URL and revoke it when the worker completes or is terminated.

#### P1 -- Important

- **P1-1** [Engineering]: No input validation on effect size — negative HR accepted (line 2518-2540)
  - For TTE endpoint, `effectToTheta` calls `Math.log(effect)`. If user enters HR=0 or HR<0, this produces -Infinity or NaN, corrupting the entire simulation.
  - Fix: Validate `effect > 0` for TTE and binary, allow any value for continuous.

- **P1-2** [Engineering]: `controlRate` input read inside `effectToTheta` but may not exist in DOM for TTE/continuous (line 2529)
  - `document.getElementById('controlRate')` returns null if the control rate input is hidden or doesn't exist for non-binary endpoints. `parseFloat(null)` = NaN.
  - Fix: Guard with `|| 0.2` fallback (already done, but should log a warning).

- **P1-3** [Statistical]: Lan-DeMets spending function formula only correct for rho >= 1 (line 1053-1057)
  - `alpha * t^rho` is the standard Lan-DeMets approximation to OBF (rho=3) and Pocock (rho=1). For rho < 1 (e.g., rho=0.5), alpha*(t) can exceed alpha at t=1, violating the spending constraint.
  - The function correctly returns alpha at t=1, but the spending increments can be negative for non-monotone spending paths.
  - Fix: Add comment that rho ∈ [1, 4] is recommended, or clamp: `return Math.min(alpha, alpha * Math.pow(t, rho))`.

- **P1-4** [Accessibility]: Tab panels lack `role="tabpanel"` and `aria-labelledby` (lines ~850-920)
  - Tabs have `role="tab"` but panels don't have `role="tabpanel"`. Screen readers can't associate panels with their tabs.
  - Fix: Add `role="tabpanel"` and `aria-labelledby="tab-{name}"` to each panel div.

- **P1-5** [Engineering]: OC summary table innerHTML builds HTML from simulation numbers (lines 2569-2585)
  - While all values are numeric (from simulation), the `scenarioLabel` is derived from user input (parsed from text field). If someone enters `<img onerror=...>` as a scenario label, it could execute.
  - Fix: Apply `escapeHtml()` to scenario labels in the OC table.

#### P2 -- Minor

- **P2-1** [Style]: Missing `<meta name="description">` and `<title>` could be more descriptive
- **P2-2** [Engineering]: Magic number `-10` used to check for "no futility" throughout (e.g., line 1207, 1297, 2466). Should be a named constant.
- **P2-3** [Statistical]: Gauss-Legendre 32-point may be insufficient for K=5 with nested recursive density. Consider adaptive quadrature or more points for K>3.
- **P2-4** [A11y]: SVG charts lack `<title>` and `<desc>` elements for screen readers.

#### False Positive Watch
- Div balance: 113/113 — CORRECT
- No Math.random() — VERIFIED (0 occurrences)
- No literal `</script>` inside JS — VERIFIED (only 1 occurrence = the actual closing tag)
- normalCDF Abramowitz & Stegun 26.2.17 — CORRECT formula
- xoshiro128** scrambler — CORRECT
- Gauss-Legendre 32-point nodes/weights — VERIFIED against published tables
