# AdaptSim

A browser-based simulator for adaptive group-sequential trial design. Reproduces rpact boundary calculations to four decimal places using a single-file HTML/JS app — no R, no Python, no install.

**Live dashboard:** <https://mahmood726-cyber.github.io/adaptsim/>

## What it does

- Group-sequential boundary calculation via Armitage-McPherson-Rowe recursive integration with 32-point Gauss-Legendre quadrature.
- Four alpha-spending functions (O'Brien-Fleming, Pocock, Hwang-Shih-DeCani, Lan-DeMets approximation).
- Binding *and* non-binding futility rules.
- Operating-characteristic Monte Carlo simulation (up to 100,000 trials).
- Protocol-ready outputs (boundary tables, ASN tables, OC plots).
- Three pre-loaded cardiovascular trial examples (DAPA-HF, EMPEROR-Reduced, SPRINT).

## Run

Open `adaptsim.html` (or `index.html`) in any modern browser. No build step.

For local development:

```bash
python -m http.server 8000
# then open http://localhost:8000/
```

## Test

```bash
python -m pytest -q
```

Tests live under `tests/`; `pytest.ini` registers the directory as a package. The suite includes:
- Boundary-accuracy parity tests against rpact reference values.
- Monte Carlo coverage checks with relaxed tolerance per `~/.claude/rules/advanced-stats.md` (atol=0.05 for stochastic simulations, NOT for deterministic estimator validation).

## Repo layout

| Path | Purpose |
|---|---|
| `adaptsim.html` | the simulator (main artifact) |
| `index.html` | landing page |
| `tests/` | pytest unit + parity tests |
| `figures/` | published figures for the manuscript |
| `manuscript_lancet_dh.md` | Lancet Digital Health submission manuscript |
| `cover_letter_lancet_dh.md` | submission cover letter |
| `e156-submission/` | E156 micro-paper bundle |
| `E156-PROTOCOL.md` | project metadata (E156 entry #3) |

## Validation

Boundary calculations match rpact (R) to 4 decimal places across all four spending functions. Monte Carlo type-I error stays within 0.002 of nominal alpha across the OC scenarios shipped with the tool.

## Key reference

Pallmann P, Bedding AW, Choodari-Oskooei B, et al. *Adaptive designs in clinical trials: why use them, and how to run and report them.* BMC Med 2018;16:29. `doi:10.1186/s12916-018-1017-7`

## License

See `LICENSE` (MIT).
