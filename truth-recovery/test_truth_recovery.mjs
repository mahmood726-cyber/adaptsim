/* test_truth_recovery.mjs — assertions for adaptsim truth-recovery.
 * Run: node truth-recovery/test_truth_recovery.mjs
 * Each test prints PASS/FAIL; process exits 1 if any fail.
 *
 * Claim under test: adaptsim's Lan-DeMets recursive-integration boundaries
 * CONTROL the type-I error at the nominal alpha across sequential looks, and
 * yield correct relative power (OBF > Pocock at matched drift). nSim is sized
 * so 3*SE comfortably brackets the nominal alpha. Seeds pinned.
 */

import { typeIError, powerCheck, boundariesFor, equalInfo, mcSE } from './harness.mjs';
import { normalQuantile } from './engine.mjs';

const NSIM = 100000;
const ALPHA = 0.025;
let failures = 0;
function check(name, cond, detail) {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${detail ? '  (' + detail + ')' : ''}`);
  if (!cond) failures++;
}

// 1. Single-look boundary equals z_alpha exactly (engine sanity).
{
  const z = boundariesFor('obf', ALPHA, equalInfo(1), false, null)[0];
  check('K=1 boundary == z_alpha', Math.abs(z - normalQuantile(1 - ALPHA)) < 1e-6,
        `z=${z.toFixed(6)}`);
}

// 2. OBF K=5 boundaries match published gsDesign one-sided values (~[4.88,3.36,2.68,2.29,2.03]).
{
  const z = boundariesFor('obf', ALPHA, equalInfo(5), false, null);
  const ref = [4.877, 3.357, 2.680, 2.290, 2.031];
  const ok = z.every((v, i) => Math.abs(v - ref[i]) < 0.01);
  check('OBF K=5 boundaries match reference', ok, `[${z.map(v => v.toFixed(3)).join(',')}]`);
}

// 3. Type-I CONTROL: OBF K=5 realised type-I within 3*SE of nominal alpha.
{
  const r = typeIError('obf', ALPHA, 5, false, NSIM, 782);
  const ok = Math.abs(r.rejectRate - ALPHA) <= 3 * mcSE(r.rejectRate, NSIM);
  check('OBF K=5 controls type-I at alpha', ok,
        `type-I=${r.rejectRate.toFixed(5)} vs ${ALPHA}`);
}

// 4. Type-I CONTROL: Pocock K=5 (the most aggressive early-spending boundary)
//    must NOT over-spend (this is exactly where a naive z/sqrt(t) approximation fails).
{
  const r = typeIError('pocock', ALPHA, 5, false, NSIM, 783);
  const ok = Math.abs(r.rejectRate - ALPHA) <= 3 * mcSE(r.rejectRate, NSIM);
  check('Pocock K=5 controls type-I at alpha (no over-spend)', ok,
        `type-I=${r.rejectRate.toFixed(5)} vs ${ALPHA}`);
}

// 5. POWER ordering: OBF retains more power than Pocock at the same drift.
{
  const theta = normalQuantile(1 - ALPHA) + normalQuantile(0.90);
  const obf = powerCheck('obf', ALPHA, 5, theta, false, NSIM, 9095).rejectRate;
  const poc = powerCheck('pocock', ALPHA, 5, theta, false, NSIM, 9095).rejectRate;
  check('OBF power > Pocock power at matched drift', obf > poc + 0.01,
        `OBF=${obf.toFixed(4)} Pocock=${poc.toFixed(4)}`);
}

console.log(failures === 0 ? '\nALL TESTS PASSED' : `\n${failures} TEST(S) FAILED`);
process.exit(failures === 0 ? 0 : 1);
