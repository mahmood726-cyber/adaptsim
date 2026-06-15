# Truth-Recovery Report — adaptsim

**Repo:** mahmood726-cyber/adaptsim
**Engine:** adaptsim.html (108 KB single-file simulator for adaptive group-sequential trials)
**Date:** 2026-06-15
**Method:** Pure boundary/spending functions extracted verbatim into engine.mjs; a
standalone seeded Brownian-motion sequential-monitoring DGP (dgp-sequential.mjs)
injects the known null (true drift = 0) and known alternatives, then applies
adaptsim's OWN computed boundaries and measures realised type-I error and power.

## Verdict: GENUINE ENGINE — TYPE-I CONTROL CONFIRMED (ship)

adaptsim is a real methods engine, not a teaching visualizer. Its boundaries are
computed by exact Lan-DeMets-style recursive numerical integration of the
canonical group-sequential joint distribution (computeBoundaries ->
buildRecursiveDensity -> computeCrossingProb, 32-point Gauss-Legendre quadrature),
NOT the naive z/sqrt(t) O'Brien-Fleming approximation. The alpha-spending function
itself is the correct OBF form alpha*(t) = 2 - 2*Phi(z_{alpha/2}/sqrt(t)).

## Measured results (nSim = 100,000, one-sided alpha = 0.025, equal information)

| Design | Boundaries (z) | Realised type-I | Nominal | Verdict |
|---|---|---|---|---|
| K=1 single look | 1.9600 | 0.02397 | 0.025 | exact (== z_alpha) |
| OBF K=3 | 3.710, 2.511, 1.993 | 0.02521 | 0.025 | CONTROLS |
| OBF K=5 | 4.877, 3.357, 2.680, 2.290, 2.031 | 0.02559 | 0.025 | CONTROLS |
| Pocock K=5 | 2.438, 2.427, 2.410, 2.397, 2.386 | 0.02493 | 0.025 | CONTROLS |

All realised type-I errors fall within 3x Monte-Carlo-SE (~0.0015) of the nominal
0.025. The OBF K=5 boundaries match published gsDesign one-sided values
(~[4.88, 3.36, 2.68, 2.29, 2.03]) to <0.01.

Power (drift theta = 3.242, ~90% single-look power):
- OBF K=5 -> 0.8933
- Pocock K=5 -> 0.8395

The OBF > Pocock power ordering at matched drift is the textbook tradeoff (OBF
spends little alpha early, preserving power for the final look), confirming the
boundaries behave correctly under H1 as well as H0.

## Contrast with the sibling TSA finding

A sibling TSA engine using the simple z_alpha/sqrt(t) OBF approximation over-spent
(type-I 0.071 at 10 looks vs nominal 0.05). adaptsim does NOT exhibit this defect
— it uses the exact recursive boundary and holds alpha across every tested look
configuration.

## Honest finding — performance, not correctness

buildRecursiveDensity re-integrates the joint density from scratch at every
density evaluation (no grid caching / forward-recursion memoization), making
boundary computation O(32^(K-1)). Measured: K=4 ~0.25 s, K=5 ~7 s, K>=6 exceeds
60 s and is effectively intractable in-browser. The math is correct; the
implementation does not scale. We therefore validated up to K=5 (the practical
ceiling) and cache boundaries once per design before simulating. The original
prompt's K=10 over/under-spend check could not be run because the engine cannot
produce K=10 boundaries in feasible time — a real limitation worth recording.

## Recommendation

- SHIP this validation. Type-I control is confirmed; the engine is statistically sound.
- Optimization (non-blocking): rewrite buildRecursiveDensity to carry a cached
  density grid forward (standard Lan-DeMets forward recursion) instead of
  re-integrating per evaluation. This changes O(32^(K-1)) -> O(K * grid^2) and
  makes K>=10 instant, enabling the high-look type-I check that currently cannot run.

## Files
- engine.mjs — extracted pure functions (PRNG, normal dist, 4 spending functions, boundary engine)
- dgp-sequential.mjs — seeded Brownian sequential-monitoring DGP
- harness.mjs — drives boundaries + DGP, reports type-I & power
- test_truth_recovery.mjs — 5 assertions (all PASS)
- REPORT.md — this file
