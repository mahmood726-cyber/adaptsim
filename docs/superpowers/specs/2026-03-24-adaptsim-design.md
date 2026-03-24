# AdaptSim — Adaptive Clinical Trial Design Simulator

## 1. Problem Statement

Adaptive clinical trial design — particularly group-sequential methods with interim analyses — is used in >90% of Phase 3 trials. Yet the tools to design these trials are locked behind expensive commercial software (FACTS ~$15K/yr, East ~$8K/yr, nQuery ~$5K/yr). Free alternatives exist in R (rpact, gsDesign) but require programming expertise. No browser-based, open-access tool exists, locking out academic groups, LMICs, and early-career trialists.

## 2. Goal

Build the world's first browser-based adaptive trial design simulator. A single-file HTML app that computes group-sequential boundaries, simulates operating characteristics via Monte Carlo, supports sample size re-estimation, and generates protocol-ready reports — all offline, free, with no server dependency.

## 3. Target

- **Primary:** Lancet Digital Health (tool paper)
- **Secondary:** F1000Research (methods paper)
- **Location:** `C:\AdaptSim\adaptsim.html`

## 4. Scope

### In scope
- Group-sequential designs with 1-5 interim analyses
- Alpha spending functions: O'Brien-Fleming, Pocock, Lan-DeMets (rho), Hwang-Shih-DeCani (gamma), custom
- Futility boundaries: binding and non-binding, beta-spending
- Monte Carlo operating characteristics (10K-100K trials)
- Sample size re-estimation (promising zone, conditional power)
- Time-to-event (HR), binary (OR/RR), and continuous (MD) endpoints
- R code export (rpact/gsDesign equivalent)
- Auto-generated methods text

### Out of scope (future phases)
- Bayesian adaptive designs (handled by MAPriors)
- Biomarker-adaptive / enrichment designs
- Platform / MAMS trials
- Response-adaptive randomization
- Dose-finding (CRM, BOIN)

## 5. Architecture

### 5.1 Single-file HTML app
Following the established pattern from 12 production apps. CSS variables, dark/light toggle, Plotly-free SVG charts, seeded PRNG, escapeHtml everywhere.

### 5.2 Tab Structure

**Tab 1: Design Builder**
- Trial parameters: endpoint type (TTE/binary/continuous), effect size, alpha (one/two-sided), target power, max sample size or events
- Interim schedule: number of looks (1-5), information fractions (equal or custom sliders)
- Efficacy spending: dropdown (OBF, Pocock, Lan-DeMets ρ, HSD γ, custom) with parameter sliders
- Futility spending: dropdown (none, binding beta-spending, non-binding) with separate spending function
- Allocation ratio slider (1:1, 2:1, 3:1)
- All inputs update boundaries in real-time (no "run" button)

**Tab 2: Boundaries**
- Primary chart: efficacy + futility boundaries vs information fraction (SVG)
- Scale toggle: Z-value / p-value / effect (HR, OR, or MD) scale
- Boundary table: look number, information fraction, z_eff, z_fut, nominal p_eff, nominal p_fut, cumulative alpha spent
- Overlay mode: compare 2 spending functions side-by-side (dropdown to select comparison)
- Visual: shaded rejection/continuation/futility regions

**Tab 3: Operating Characteristics**
- Scenario panel: define 3-8 true effect sizes (null + alternatives)
- Simulation control: number of simulations (10K/50K/100K), seed
- Run button → Web Worker Monte Carlo (progress bar, cancellable)
- Results dashboard:
  - Power curve: power vs true effect size
  - ASN curve: expected sample size vs true effect
  - Stopping probabilities: stacked bar at each look per scenario
  - Conditional power distribution at each interim
- Summary table: Type I error, power, ASN, median sample size, probability of early stopping per scenario

**Tab 4: Sample Size Re-estimation (SSR)**
- SSR trigger: at which interim look (dropdown)
- Method: conditional power-based (Cui-Hung-Wang) or promising zone (Chen-DeMets-Lan)
- Target conditional power for SSR (default 80%)
- Promising zone bounds: CP_low (e.g., 0.36) to CP_high (e.g., 0.80)
- Max sample size cap (fold-increase limit, e.g., 2x)
- Blinded vs unblinded toggle
- OC comparison: fixed design vs SSR-adaptive design (overlay curves)

**Tab 5: Report & Export**
- Auto-generated methods paragraph (SPIRIT/ICH E9 R1 compliant)
- Boundary table formatted for protocol appendix
- OC summary table
- R code export: equivalent rpact::getDesignGroupSequential() + rpact::getSimulationResults()
- Print-ready layout (A4 CSS)
- Copy buttons for each section

### 5.3 Core Mathematical Engine

#### Alpha spending functions

| Function | Formula | Parameter |
|----------|---------|-----------|
| O'Brien-Fleming | α*(t) = 2 - 2Φ(z_{α/2} / √t) | none |
| Pocock | α*(t) = α · ln(1 + (e-1)·t) | none |
| Lan-DeMets ρ | α*(t) = α · t^ρ | ρ ∈ [0.5, 4] (ρ=1 Pocock-like, ρ=3 OBF-like) |
| Hwang-Shih-DeCani | α*(t) = α · (1 - e^(-γt)) / (1 - e^(-γ)) | γ ∈ [-10, 10] (γ→0 = Pocock, γ=-4 ≈ OBF) |
| Custom | User-specified α*(t_k) at each look | direct input |

#### Boundary computation (recursive numerical integration)

At each look k (k = 1, ..., K):
1. Compute incremental alpha: Δα_k = α*(t_k) - α*(t_{k-1})
2. Find z_k such that P(reject at look k | survived looks 1..k-1, H0) = Δα_k
3. This requires computing the joint probability under the canonical joint distribution of sequential Z-statistics:
   - Z = (Z_1, ..., Z_K) ~ MVN(0, Σ) under H0
   - Σ_{ij} = √(t_i / t_j) for i ≤ j (the independent increments covariance)
4. Use recursive numerical integration (Armitage-McPherson-Rowe 1969):
   - At look 1: z_1 = Φ^{-1}(1 - Δα_1/2) for two-sided
   - At look k>1: integrate over the continuation region at looks 1..k-1, find z_k by bisection

Implementation: numerical integration via Gauss-Legendre quadrature (32-64 points). This is exact to ~6 decimal places and fast enough for real-time updates.

#### Monte Carlo simulation engine

For each simulated trial under true effect δ:
1. Generate sequential test statistics: Z_k ~ N(δ√(t_k · I_max), 1) with independent increments
2. At each look, compare Z_k to efficacy boundary (reject) and futility boundary (stop for futility)
3. Record: which look stopped, final Z, final decision

Use Web Workers to run simulations off the main thread. Batch in chunks of 1,000 for responsive progress reporting.

#### Sample size re-estimation

At interim look m with observed Z_m:
1. Compute conditional power: CP = 1 - Φ((z_final - Z_m·√(t_m)) / √(1-t_m) - δ·√(I_remaining))
2. If CP ∈ [CP_low, CP_high] (promising zone): increase n to achieve target CP
3. New n = n_original · (z_{1-β} + z_{α/2})² / (Z_m·√(t_m) + z_{1-β}·√(1-t_m))²
4. Cap at max_fold · n_original
5. Adjust final boundary using Cui-Hung-Wang (CHW) approach to preserve Type I error

#### Endpoint-specific conversions

| Endpoint | Test statistic | Effect → Z |
|----------|---------------|------------|
| Time-to-event | log-rank Z | Z = log(HR) · √(d/4) where d = events |
| Binary | Z-test on proportions | Z = (p1-p0) / √(p̄(1-p̄)(1/n1+1/n2)) |
| Continuous | t-test on means | Z = (μ1-μ0) / √(σ²(1/n1+1/n2)) |

Information fraction for TTE: t_k = d_k / d_max (events observed / total events planned).

### 5.4 Built-in Examples

1. **DAPA-HF style** — TTE, HR=0.74, 2 interim looks, OBF spending, non-binding futility
2. **EMPEROR-Reduced style** — TTE, HR=0.75, 1 interim look, Lan-DeMets ρ=3
3. **SPRINT style** — Binary (composite), OR=0.75, 3 interim looks with early stopping for efficacy achieved at look 2

## 6. Validation Strategy

1. **rpact cross-validation**: For each built-in example, generate equivalent rpact R code and compare boundaries, power, and ASN to 4 decimal places
2. **Known boundary tables**: Verify against published Lan-DeMets tables (Jennison & Turnbull 2000, Table 2.4)
3. **Type I error verification**: Under H0, simulated rejection rate must be ≤ α + 0.002 (simulation error tolerance)
4. **Selenium tests**: 30+ tests covering all tabs, edge cases (k=1, extreme gamma, zero futility)

## 7. Non-Goals

- NOT a sample size calculator (focused on adaptive design operating characteristics)
- NOT for Bayesian adaptive designs (MAPriors covers that)
- NOT for dose-finding or biomarker-adaptive designs
- NOT for regulatory submission document generation (only methods text)

## 8. Success Criteria

1. Boundaries match rpact within ±0.001 on all 3 built-in examples
2. Monte Carlo Type I error ≤ α + 0.002 with 100K simulations
3. Real-time boundary updates (<100ms response to slider changes)
4. Web Worker simulation completes 100K trials in <10 seconds
5. All 30+ Selenium tests pass
6. Methods text passes ICH E9 R1 checklist
