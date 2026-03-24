# AdaptSim Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the world's first browser-based adaptive trial design simulator as a single-file HTML app.

**Architecture:** Single HTML file with embedded CSS + JS. Core math engine (spending functions, boundary computation via Gauss-Legendre quadrature, Monte Carlo via Web Workers). 5 tabs: Design, Boundaries, OC, SSR, Report. SVG charts, no external dependencies.

**Tech Stack:** HTML5, CSS3 (variables, dark mode), vanilla JavaScript, SVG, Web Workers (inline blob), scipy-equivalent math in JS.

---

### Task 1: HTML Shell + Tab Navigation + CSS Theme

**Files:**
- Create: `C:\AdaptSim\adaptsim.html`

- [ ] Build the full HTML skeleton: 5 tabs, CSS variables, dark/light toggle, responsive layout
- [ ] Tab 1 (Design): input form with endpoint type, effect size, alpha, power, K looks, info fractions, spending function dropdowns, allocation ratio
- [ ] Tab 2 (Boundaries): empty SVG container + boundary table placeholder
- [ ] Tab 3 (OC): scenario inputs + simulation controls + results containers
- [ ] Tab 4 (SSR): SSR method, trigger look, promising zone inputs
- [ ] Tab 5 (Report): methods text + export containers
- [ ] Tab switching JS + 3 built-in example buttons
- [ ] Commit

### Task 2: Math Engine — Spending Functions + Normal Distribution

- [ ] Implement in JS: normalCDF, normalQuantile (Beasley-Springer-Moro), tCDF, tQuantile (via scipy-equivalent)
- [ ] Implement 4 spending functions: OBF, Pocock, Lan-DeMets(rho), HSD(gamma)
- [ ] Each returns alpha*(t) given total alpha and information fraction t
- [ ] Commit

### Task 3: Boundary Computation Engine

- [ ] Implement Gauss-Legendre quadrature (32 nodes + weights)
- [ ] Implement recursive boundary search: at look k, find z_k via bisection such that cumulative rejection probability = alpha*(t_k)
- [ ] The key algorithm: integrate over continuation region at looks 1..k-1 using the independent increments property of sequential Z-statistics
- [ ] Support two-sided tests (symmetric boundaries)
- [ ] Support futility boundaries (beta-spending analog)
- [ ] Return: array of {look, info_frac, z_eff, z_fut, nominal_p_eff, cum_alpha_spent}
- [ ] Commit

### Task 4: Boundary Visualization (Tab 2)

- [ ] SVG chart: x-axis = information fraction (0-1), y-axis = Z-value
- [ ] Plot efficacy boundary (upper), futility boundary (lower if enabled)
- [ ] Shaded regions: rejection (above efficacy), futility (below futility), continuation (between)
- [ ] Scale toggle: Z / p-value / effect scale
- [ ] Boundary table below chart
- [ ] Overlay mode: compare 2 spending functions
- [ ] Wire up real-time updates from Tab 1 sliders
- [ ] Commit

### Task 5: Monte Carlo OC Engine (Web Worker)

- [ ] Inline Web Worker (Blob URL) that simulates group-sequential trials
- [ ] Input: boundaries, info fractions, true effect (delta), n_sims, seed
- [ ] For each trial: generate sequential Z-stats with drift delta*sqrt(t_k*I_max), check boundaries at each look
- [ ] Output: rejection count, futility stop count, look-specific stop counts, total sample sizes
- [ ] Progress messages back to main thread every 1000 sims
- [ ] Cancellation support via Worker.terminate()
- [ ] Commit

### Task 6: OC Dashboard (Tab 3)

- [ ] Scenario panel: add/remove scenarios (null + up to 7 alternatives)
- [ ] Run button → dispatch to Web Worker → progress bar
- [ ] Power curve SVG: power vs true effect
- [ ] ASN curve SVG: expected sample size vs true effect
- [ ] Stopping probability stacked bars at each look
- [ ] Summary table: Type I error, power, ASN, median N, P(early stop) per scenario
- [ ] Commit

### Task 7: Sample Size Re-estimation (Tab 4)

- [ ] Conditional power formula implementation
- [ ] Promising zone logic (Chen-DeMets-Lan)
- [ ] SSR calculation: new n given observed Z at interim
- [ ] Blinded SSR approximation
- [ ] OC comparison: overlay fixed vs SSR-adaptive on same chart
- [ ] Commit

### Task 8: Report & Export (Tab 5)

- [ ] Auto-generated methods text from current design parameters
- [ ] Boundary table formatted for copy-paste into protocol
- [ ] R code export: rpact::getDesignGroupSequential() equivalent
- [ ] Print-ready CSS (@media print)
- [ ] Copy-to-clipboard buttons
- [ ] Commit

### Task 9: Built-in Examples + Validation

- [ ] DAPA-HF: TTE HR=0.74, K=2, OBF, non-binding futility
- [ ] EMPEROR-Reduced: TTE HR=0.75, K=1, Lan-DeMets rho=3
- [ ] SPRINT: binary composite OR=0.75, K=3
- [ ] Verify boundaries against rpact reference values
- [ ] Commit

### Task 10: Testing + Final Polish

- [ ] Selenium test suite: 30+ tests across all tabs
- [ ] Div balance check
- [ ] Accessibility: keyboard nav, ARIA labels, focus management
- [ ] Dark mode verification
- [ ] Final commit
