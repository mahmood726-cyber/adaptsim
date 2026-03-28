Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

AdaptSim: A Browser-Based Simulator for Adaptive Group-Sequential Trial Design Matching rpact to Four Decimal Places

Can a browser-based tool reproduce adaptive group-sequential trial designs with the accuracy of established software while remaining usable at the point of care? We implemented the Armitage-McPherson-Rowe recursive integration algorithm with 32-point Gauss-Legendre quadrature in a single-file browser application supporting four alpha-spending functions and both binding and non-binding futility rules. AdaptSim calculates monitoring boundaries, simulates operating characteristics with up to 100,000 Monte Carlo trials, and exports protocol-ready outputs with matching R code. Boundary calculations matched rpact to four decimal places, with mean sensitivity to true effect of 0.98 (95% CI 0.96 to 0.99) and Monte Carlo type I error within 0.002 of nominal alpha. Three worked examples based on DAPA-HF, EMPEROR-Reduced, and SPRINT showed that the interface can support realistic cardiovascular trial scenarios without local installation. A free browser implementation could widen access to adaptive design methods that otherwise depend on specialist software and scripting skills. The current version is limited to group-sequential designs and does not yet support adaptive enrichment or platform trials.

Outside Notes

Type: methods
Primary estimand: Boundary computation accuracy vs rpact (decimal places of agreement)
App: AdaptSim v1.0
Data: Validated against rpact R package; 3 pre-loaded CV trial examples
Code: https://github.com/mahmood726-cyber/adaptsim
Version: 1.0
Certainty: high
Validation: Author reviewed draft

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
