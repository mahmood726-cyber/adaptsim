[DATE]

The Editor
The Lancet Digital Health

Dear Editor,

We submit the manuscript **"AdaptSim: A Zero-Install Browser-Based Adaptive Trial Design Simulator with Group-Sequential Boundaries and Sample Size Re-Estimation"** for consideration in *The Lancet Digital Health*.

**The gap:** Adaptive trial design — group-sequential monitoring, alpha-spending boundaries, conditional power-based sample size re-estimation — requires commercial software (EAST, nQuery) or R programming expertise (rpact, gsDesign). This creates a barrier for clinical investigators who need to explore design options during protocol development.

**Our solution:** AdaptSim is the first fully browser-based tool for adaptive trial design. It computes O'Brien-Fleming, Pocock, Lan-DeMets (alpha-spending), and Hwang-Shih-DeCani boundaries using 32-point Gauss-Legendre quadrature, runs Monte Carlo operating characteristics with a seeded PRNG for reproducibility, and exports equivalent R code for regulatory documentation.

**Validation:** Boundary computations match the rpact R package to 4 decimal places across all four alpha-spending functions. Type I error rates in 100,000-simulation Monte Carlo runs are within 0.002 of nominal alpha. Three pre-loaded examples from landmark trials (DAPA-HF, EMPEROR-Reduced, SPRINT) demonstrate real-world applicability.

**Lancet Digital Health fit:** This tool democratises access to methodology that currently requires expensive commercial software or statistical programming expertise, directly enabling better-designed clinical trials.

The manuscript has not been submitted elsewhere. The tool and code are freely available.

Yours sincerely,

Mahmood Ahmad
Royal Free Hospital, London, UK
mahmood.ahmad2@nhs.net
