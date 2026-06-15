/* harness.mjs — truth-recovery driver for adaptsim.
 *
 * Q1 (type-I control): Under the null (theta=0), do adaptsim's OWN computed
 *     OBF / Pocock boundaries hold the cumulative type-I error at alpha across
 *     all planned looks? The sibling TSA's naive z/sqrt(t) OBF over-spent
 *     (type-I 0.071 at 10 looks). adaptsim instead uses the exact Lan-DeMets
 *     recursive-integration boundary, which should control type-I at alpha.
 *
 * Q2 (power): Under a fixed H1 drift, realised power is reported (sanity).
 *
 * PERFORMANCE FINDING: adaptsim's buildRecursiveDensity re-integrates from
 *     scratch at every density evaluation (no grid caching), so boundary
 *     computation is O(32^(K-1)). K=5 ~7s; K>=6 is intractable (>60s). We
 *     therefore test up to K=5 and document this limit. Boundaries are computed
 *     ONCE per design and cached, then reused across all Monte-Carlo replicates.
 */

import { computeBoundaries, getSpendingFunction, normalQuantile } from './engine.mjs';
import { simulate } from './dgp-sequential.mjs';

export function equalInfo(K) {
  const a = [];
  for (let k = 1; k <= K; k++) a.push(k / K);
  return a;
}

const _bcache = new Map();
export function boundariesFor(spendType, alpha, infoFracs, twoSided, param) {
  const key = `${spendType}|${alpha}|${infoFracs.join(',')}|${twoSided}|${param}`;
  if (_bcache.has(key)) return _bcache.get(key);
  const spendFn = getSpendingFunction(spendType, param);
  const b = computeBoundaries(alpha, infoFracs, spendFn, twoSided, 'none', null, 0.2);
  const zEff = b.map(x => x.z_eff);
  _bcache.set(key, zEff);
  return zEff;
}

export function typeIError(spendType, alpha, K, twoSided, nSim, seed, param) {
  const infoFracs = equalInfo(K);
  const zEff = boundariesFor(spendType, alpha, infoFracs, twoSided, param);
  const res = simulate(infoFracs, zEff, 0.0, twoSided, nSim, seed);
  return { zEff, infoFracs, ...res };
}

export function powerCheck(spendType, alpha, K, theta, twoSided, nSim, seed, param) {
  const infoFracs = equalInfo(K);
  const zEff = boundariesFor(spendType, alpha, infoFracs, twoSided, param);
  const res = simulate(infoFracs, zEff, theta, twoSided, nSim, seed);
  return { zEff, theta, ...res };
}

export function mcSE(p, n) { return Math.sqrt(p * (1 - p) / n); }

function main() {
  const NSIM = 100000;
  const ALPHA = 0.025;       // one-sided default in adaptsim UI
  const lines = [];
  lines.push('=== adaptsim truth-recovery: type-I & power ===');
  lines.push(`nSim=${NSIM}, alpha=${ALPHA} (one-sided), equal info fractions`);
  lines.push('');

  // Single-look sanity
  {
    const r = typeIError('obf', ALPHA, 1, false, NSIM, 12345);
    lines.push(`[K=1] zEff=${r.zEff[0].toFixed(4)} (expect ${normalQuantile(1-ALPHA).toFixed(4)}), type-I=${r.rejectRate.toFixed(5)} +/- ${(3*mcSE(r.rejectRate,NSIM)).toFixed(5)}`);
    lines.push('');
  }

  for (const [name, K] of [['obf',3],['obf',5],['pocock',5]]) {
    const r = typeIError(name, ALPHA, K, false, NSIM, 777 + K);
    const se = mcSE(r.rejectRate, NSIM);
    const verdict = Math.abs(r.rejectRate - ALPHA) <= 3 * se ? 'CONTROLS' : (r.rejectRate > ALPHA ? 'OVER-SPENDS' : 'UNDER-SPENDS');
    lines.push(`[${name.toUpperCase()} K=${K}] boundaries=[${r.zEff.map(z=>z.toFixed(3)).join(', ')}]`);
    lines.push(`            type-I=${r.rejectRate.toFixed(5)} (nominal ${ALPHA}, 3*SE=${(3*se).toFixed(5)}) -> ${verdict}`);
  }
  lines.push('');

  const thetaTarget = normalQuantile(1 - ALPHA) + normalQuantile(0.90); // ~90% single-look power
  for (const [name, K] of [['obf',5],['pocock',5]]) {
    const r = powerCheck(name, ALPHA, K, thetaTarget, false, NSIM, 9090 + K);
    lines.push(`[${name.toUpperCase()} K=${K}] power at theta=${thetaTarget.toFixed(3)} -> ${r.rejectRate.toFixed(4)}`);
  }

  console.log(lines.join('\n'));
}

if (process.argv[1] && process.argv[1].endsWith('harness.mjs')) {
  main();
}
