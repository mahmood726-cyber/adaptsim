/* dgp-sequential.mjs — STANDALONE seeded known-truth DGP for adaptsim.
 *
 * Generates the canonical group-sequential test-statistic process under a
 * specified true drift, then applies adaptsim's OWN computed efficacy boundaries
 * and measures the realised type-I error (drift=0) / power (drift>0) across the
 * planned looks.
 *
 * Canonical structure (Jennison & Turnbull / Lan-DeMets):
 *   Brownian motion B(t) with drift theta and independent increments.
 *   At information fraction t_k, the score is B(t_k) = theta*t_k + W(t_k),
 *   W standard Brownian motion. The standardized statistic is
 *     Z_k = B(t_k) / sqrt(t_k)
 *   so under H0 (theta=0): Z_k ~ N(0, 1), Cov(Z_i,Z_j) = sqrt(t_i/t_j), i<=j.
 *   Under H1: E[Z_k] = theta * sqrt(t_k)  (theta in score-scale "drift").
 *
 * We simulate by drawing independent increments dW_j ~ N(0, t_j - t_{j-1}),
 * accumulating B, and standardizing — this reproduces the exact correlation
 * structure that adaptsim's boundaries are derived to control.
 */

import { seedPRNG, normalRandom } from './engine.mjs';

/**
 * Run N replicates of the sequential process.
 * @param {number[]} infoFracs  increasing information fractions in (0,1], last = 1
 * @param {number[]} zEff       efficacy boundaries (one per look) from adaptsim
 * @param {number}   theta      true drift (0 = null). Power scale: E[Z_K]=theta.
 *                              (theta is the standardized effect at full information)
 * @param {boolean}  twoSided   reject if |Z|>zEff (true) or Z>zEff (false)
 * @param {number}   nSim
 * @param {number}   seed
 * @returns {{rejectRate:number, rejectByLook:number[], firstLookRate:number[], nSim:number}}
 */
export function simulate(infoFracs, zEff, theta, twoSided, nSim, seed) {
  const K = infoFracs.length;
  const rng = seedPRNG(seed >>> 0);
  // increments of the standard Brownian motion W over each segment
  const dt = infoFracs.map((t, k) => k === 0 ? t : t - infoFracs[k - 1]);

  let rejected = 0;
  const rejectByLook = new Array(K).fill(0); // count of rejections AT each look (first crossing)

  for (let s = 0; s < nSim; s++) {
    let B = 0;            // accumulated Brownian path WITH drift, in score scale
    let crossed = false;
    for (let k = 0; k < K; k++) {
      // increment of W over [t_{k-1}, t_k] has variance dt[k]
      const dW = normalRandom(rng) * Math.sqrt(dt[k]);
      // drift contribution: theta is standardized effect at t=1 (E[Z_K]=theta),
      // score-scale drift per unit info = theta, so dB_drift = theta*dt
      B += theta * dt[k] + dW;
      const tk = infoFracs[k];
      const Z = B / Math.sqrt(tk);
      const hit = twoSided ? (Math.abs(Z) > zEff[k]) : (Z > zEff[k]);
      if (hit) {
        rejected++;
        rejectByLook[k]++;
        crossed = true;
        break; // sequential: stop at first crossing
      }
    }
    void crossed;
  }

  return {
    rejectRate: rejected / nSim,
    rejectByLook: rejectByLook.map(c => c / nSim),
    nSim,
  };
}
