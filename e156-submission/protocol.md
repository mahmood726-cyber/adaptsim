Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Protocol: AdaptSim: A Browser-Based Simulator for Adaptive Group-Sequential Trial Design Matching rpact to Four Decimal Places

This protocol describes a methods-validation study of AdaptSim, a browser-based tool for adaptive group-sequential trial design. We will compare boundary calculations from the application with outputs from the rpact R package across predefined design scenarios spanning four alpha-spending functions and both binding and non-binding futility rules. Primary outcome is decimal-place agreement in efficacy and futility boundaries, with secondary outcomes including type I error, power, and run time across Monte Carlo simulations of up to 100,000 trials. Validation scenarios will include worked examples modelled on DAPA-HF, EMPEROR-Reduced, and SPRINT to assess practical performance in cardiovascular trial settings. Analyses will be deterministic where possible, with simulation seeds fixed and all code versioned for reproducibility. Results will be reported with confidence intervals and accompanied by matching R code for external verification. The study is limited to group-sequential designs and will not evaluate adaptive enrichment, response-adaptive randomisation, or platform-trial architectures.

Outside Notes

Type: protocol
Primary estimand: Boundary computation accuracy vs rpact (decimal places of agreement)
App: AdaptSim v1.0
Code: https://github.com/mahmood726-cyber/adaptsim
Date: 2026-03-26
Validation: Author reviewed draft

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

LLM assistance was used for drafting and language editing. The author reviewed and edited the manuscript and takes responsibility for the final content.



