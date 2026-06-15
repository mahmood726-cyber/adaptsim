/* engine.mjs — pure functions extracted VERBATIM from adaptsim.html (truth-recovery sweep).
   Source: adaptsim.html lines ~935-1380. ADDITIVE; no edits to original logic. */

function xoshiro128ss(a, b, c, d) {
  return function() {
    var t = b << 9, r = b * 5;
    r = (r << 7 | r >>> 25) * 9;
    c ^= a; d ^= b; b ^= c; a ^= d;
    c ^= t;
    d = d << 11 | d >>> 21;
    return (r >>> 0) / 4294967296;
  };
}

function seedPRNG(seed) {
  /* splitmix32 to generate initial state from a single seed */
  function sm32(a) {
    return function() {
      a |= 0; a = a + 0x9e3779b9 | 0;
      var t = a ^ a >>> 16; t = Math.imul(t, 0x21f0aaad);
      t = t ^ t >>> 15; t = Math.imul(t, 0x735a2d97);
      return ((t = t ^ t >>> 15) >>> 0);
    };
  }
  var sm = sm32(seed);
  return xoshiro128ss(sm(), sm(), sm(), sm());
}

function clamp(v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; }
function fmt(v, d) { return v == null ? '-' : v.toFixed(d != null ? d : 4); }

/* ===== 2. NORMAL DISTRIBUTION ===== */

/* Standard normal CDF via rational approximation (Abramowitz & Stegun 26.2.17) */
function normalCDF(x) {
  if (x < -8) return 0;
  if (x > 8) return 1;
  var neg = x < 0;
  if (neg) x = -x;
  var p = 0.2316419;
  var b1 = 0.319381530, b2 = -0.356563782, b3 = 1.781477937, b4 = -1.821255978, b5 = 1.330274429;
  var t = 1.0 / (1.0 + p * x);
  var t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
  var pdf = Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
  var cdf = 1.0 - pdf * (b1 * t + b2 * t2 + b3 * t3 + b4 * t4 + b5 * t5);
  return neg ? 1 - cdf : cdf;
}

/* Standard normal PDF */
function normalPDF(x) {
  return Math.exp(-0.5 * x * x) / Math.sqrt(2 * Math.PI);
}

/* Inverse normal (Beasley-Springer-Moro algorithm) */
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  var a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  var b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  var c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
     4.374664141464968e+00, 2.938163982698783e+00
  ];
  var d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];

  var pLow = 0.02425, pHigh = 1 - pLow;
  var q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5]) * q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

/* Box-Muller normal random using seeded PRNG */
function normalRandom(rng) {
  var u1 = rng(), u2 = rng();
  while (u1 === 0) u1 = rng();
  return Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
}


/* ===== 3. ALPHA SPENDING FUNCTIONS ===== */

function alphaSpendOBF(t, alpha) {
  if (t <= 0) return 0;
  if (t >= 1) return alpha;
  var za2 = normalQuantile(1 - alpha / 2);
  return 2 - 2 * normalCDF(za2 / Math.sqrt(t));
}

function alphaSpendPocock(t, alpha) {
  if (t <= 0) return 0;
  if (t >= 1) return alpha;
  return alpha * Math.log(1 + (Math.E - 1) * t);
}

function alphaSpendLanDeMets(t, alpha, rho) {
  if (t <= 0) return 0;
  if (t >= 1) return alpha;
  return Math.min(alpha, alpha * Math.pow(t, rho));
}

function alphaSpendHSD(t, alpha, gamma) {
  if (t <= 0) return 0;
  if (t >= 1) return alpha;
  if (Math.abs(gamma) < 1e-8) {
    return alpha * t;
  }
  return alpha * (1 - Math.exp(-gamma * t)) / (1 - Math.exp(-gamma));
}

function getSpendingFunction(type, param) {
  switch (type) {
    case 'obf': return function(t, a) { return alphaSpendOBF(t, a); };
    case 'pocock': return function(t, a) { return alphaSpendPocock(t, a); };
    case 'landemets': return function(t, a) { return alphaSpendLanDeMets(t, a, param); };
    case 'hsd': return function(t, a) { return alphaSpendHSD(t, a, param); };
    default: return function(t, a) { return alphaSpendOBF(t, a); };
  }
}

/* Gauss-Legendre quadrature nodes and weights (32-point) */
var GL32_NODES = [], GL32_WEIGHTS = [];
(function() {
  /* Pre-computed 32-point Gauss-Legendre on [-1,1] */
  var n = [
    0.0483076656877383,0.1444719615827965,0.2392873622521371,0.3318686022821276,
    0.4213512761306353,0.5068999089322294,0.5877157572407623,0.6630442669302152,
    0.7321821187402897,0.7944837959679424,0.8493676137325700,0.8963211557660521,
    0.9349060759377397,0.9647622555875064,0.9856115115452684,0.9972638618494816
  ];
  var w = [
    0.0965400885147278,0.0956387200792749,0.0938443990808046,0.0911738786957639,
    0.0876520930044038,0.0833119242269467,0.0781938957870703,0.0723457941088485,
    0.0658222227763618,0.0586840934785355,0.0509980592623762,0.0428358980222267,
    0.0342738629130214,0.0253920653092621,0.0162743947309057,0.0070186100094701
  ];
  for (var i = n.length - 1; i >= 0; i--) {
    GL32_NODES.push(-n[i]);
    GL32_WEIGHTS.push(w[i]);
  }
  for (var i = 0; i < n.length; i++) {
    GL32_NODES.push(n[i]);
    GL32_WEIGHTS.push(w[i]);
  }
})();

/**
 * Numerical integration using 32-point Gauss-Legendre quadrature.
 * Integrates fn(x) from a to b.
 */
function gaussLegendre(fn, a, b) {
  var mid = (a + b) / 2;
  var half = (b - a) / 2;
  var sum = 0;
  for (var i = 0; i < 32; i++) {
    sum += GL32_WEIGHTS[i] * fn(mid + half * GL32_NODES[i]);
  }
  return sum * half;
}

/**
 * Compute group-sequential boundaries using recursive numerical integration.
 *
 * The canonical joint distribution: Z_k = sum of independent increments.
 * Under H0, increments W_k ~ N(0, t_k - t_{k-1}).
 * Z_k = (sum_{j=1}^{k} W_j) / sqrt(t_k)
 * Cov(Z_i, Z_j) = sqrt(t_i/t_j) for i <= j (independent increments structure).
 *
 * The algorithm proceeds recursively: at look k we need
 *   P(Z_k > c_k, Z_1 < c_1, ..., Z_{k-1} < c_{k-1}) = delta_alpha_k
 *
 * Using conditional distributions:
 *   Z_k | Z_{k-1}=z ~ N(z * sqrt(t_{k-1}/t_k), 1 - t_{k-1}/t_k)
 *
 * We carry forward the density of Z_{k-1} on the continuation region and integrate.
 */
function computeBoundaries(alpha, infoFracs, spendFn, twoSided, futType, futSpendFn, beta) {
  var K = infoFracs.length;
  var boundaries = [];
  var alphaSpentPrev = 0;
  var betaSpentPrev = 0;
  var hasFutility = futType !== 'none';
  var effectiveBeta = beta != null ? beta : 0.2;

  /* At each look we store {z_eff, z_fut, nomP_eff, nomP_fut, cumAlpha, cumBeta} */

  /* The density array: we carry forward the joint density evaluated on a grid */
  /* prevDensity[i] = density at grid point z_i in continuation region at look k-1 */

  for (var k = 0; k < K; k++) {
    var t = infoFracs[k];
    var targetAlpha = spendFn(t, alpha);
    var deltaAlpha = targetAlpha - alphaSpentPrev;
    if (deltaAlpha < 1e-15) deltaAlpha = 1e-15;

    var targetBeta = 0, deltaBeta = 0;
    if (hasFutility && k < K - 1) {
      targetBeta = futSpendFn(t, effectiveBeta);
      deltaBeta = targetBeta - betaSpentPrev;
      if (deltaBeta < 1e-15) deltaBeta = 1e-15;
    }

    var zEff, zFut;

    if (k === 0) {
      /* First look: simple quantile */
      if (twoSided) {
        zEff = normalQuantile(1 - deltaAlpha / 2);
      } else {
        zEff = normalQuantile(1 - deltaAlpha);
      }

      /* Futility at first look */
      if (hasFutility && k < K - 1) {
        /* For non-binding futility, find z_fut such that under H0,
           P(Z_1 < z_fut) = deltaBeta (using beta spending on the lower tail) */
        /* Under H1: P(Z < z_fut | theta) should be deltaBeta
           But for simplicity and standard practice we define futility z under H0 */
        /* Standard approach: use the beta spending to define futility */
        zFut = normalQuantile(deltaBeta);
      } else {
        zFut = -Infinity;
      }
    } else {
      var tPrev = infoFracs[k - 1];
      var rho = Math.sqrt(tPrev / t);
      var sigma = Math.sqrt(1 - tPrev / t);

      /* Build the continuation region density from previous look */
      var prevBounds = boundaries[k - 1];
      /* For two-sided tests, the lower boundary of the continuation region
         is -z_eff (symmetric rejection), unless an explicit futility boundary
         is set that is higher than -z_eff */
      var prevLo;
      if (twoSided) {
        prevLo = -prevBounds.z_eff;
        if (prevBounds.z_fut > -10 && prevBounds.z_fut > prevLo) {
          prevLo = prevBounds.z_fut;
        }
      } else {
        prevLo = prevBounds.z_fut > -10 ? prevBounds.z_fut : -8;
      }
      var prevHi = prevBounds.z_eff;

      /* For k >= 2, we need the joint probability by integrating over the continuation region.
         We use a recursive approach: carry the density of the test statistic
         in the continuation region forward.

         At look k, the crossing probability:
         P(Z_k > c_k, survived) = integral over z in [lo, hi] of
           g(z) * (1 - Phi((c_k - rho*z) / sigma)) dz

         where g(z) is the density of Z_{k-1} in the continuation region,
         accounting for all previous boundaries.

         For look 2 (k=1), g(z) = phi(z) (standard normal) restricted to continuation region.
         For look k > 2, g(z) is the recursively computed density from look k-1.
      */

      /* We need to compute the recursive density. For simplicity and accuracy,
         we carry the density as a function that recursively integrates. */

      /* Actually, for a practical implementation with up to 5 looks,
         we can use nested quadrature. For look k, we need (k-1)-dimensional
         integration. With K<=5, this is at most 4-dimensional, which is
         feasible with moderate quadrature points (fewer for inner dimensions). */

      /* However, the standard efficient approach is the recursive density method:
         At look j, define:
           f_1(z) = phi(z) for z in [z_fut_1, z_eff_1]
           f_j(z) = integral over z' in [z_fut_{j-1}, z_eff_{j-1}] of
                     f_{j-1}(z') * phi_conditional(z | z') dz'
         where phi_conditional(z | z') = N(rho_j * z', sigma_j^2)

         Then P(cross at look k) = integral of f_{k-1}(z) * P(Z_k > c_k | z) dz
      */

      /* Build recursive density function */
      var densityFn = buildRecursiveDensity(k, boundaries, infoFracs, twoSided);

      /* Bisection to find z_eff */
      var zLo = 0, zHi = 8;

      /* Binary search for the z_eff such that crossing prob = deltaAlpha */
      for (var iter = 0; iter < 80; iter++) {
        var zMid = (zLo + zHi) / 2;
        var prob = computeCrossingProb(zMid, densityFn, prevLo, prevHi, rho, sigma, twoSided);
        if (prob > deltaAlpha) {
          zLo = zMid;
        } else {
          zHi = zMid;
        }
      }
      zEff = (zLo + zHi) / 2;

      /* Futility boundary for this look */
      if (hasFutility && k < K - 1) {
        /* Find z_fut such that under H0, P(Z_k < z_fut, survived) = deltaBeta */
        /* This uses the same recursive density but integrating the lower tail */
        var fLo = -8, fHi = zEff - 0.01;
        for (var iter = 0; iter < 80; iter++) {
          var fMid = (fLo + fHi) / 2;
          var fProb = computeFutilityCrossingProb(fMid, densityFn, prevLo, prevHi, rho, sigma);
          if (fProb < deltaBeta) {
            fLo = fMid;
          } else {
            fHi = fMid;
          }
        }
        zFut = (fLo + fHi) / 2;
      } else if (k === K - 1 && hasFutility) {
        /* At final look, futility = efficacy (reject or fail to reject) */
        zFut = zEff;
      } else {
        zFut = -Infinity;
      }
    }

    var nomPEff = twoSided ? 2 * (1 - normalCDF(zEff)) : 1 - normalCDF(zEff);

    boundaries.push({
      look: k + 1,
      infoFrac: t,
      z_eff: zEff,
      z_fut: zFut,
      nomP_eff: nomPEff,
      nomP_fut: zFut > -10 ? normalCDF(zFut) : 0,
      cumAlpha: targetAlpha,
      cumBeta: targetBeta
    });

    alphaSpentPrev = targetAlpha;
    if (hasFutility) betaSpentPrev = targetBeta;
  }

  return boundaries;
}

/**
 * Build recursive density function for look k.
 * Returns a function f(z) representing the density of Z_{k-1} in the continuation region.
 */
function buildRecursiveDensity(k, boundaries, infoFracs, twoSided) {
  if (k === 1) {
    /* Density at look 1 is just the standard normal, restricted to continuation region */
    return function(z) {
      return normalPDF(z);
    };
  }

  /* For k > 1, we recursively integrate */
  var prevDensity = buildRecursiveDensity(k - 1, boundaries, infoFracs, twoSided);
  var prevBounds = boundaries[k - 2]; /* k-2 because boundaries is 0-indexed */
  var prevLo;
  if (twoSided) {
    prevLo = -prevBounds.z_eff;
    if (prevBounds.z_fut > -10 && prevBounds.z_fut > prevLo) {
      prevLo = prevBounds.z_fut;
    }
  } else {
    prevLo = prevBounds.z_fut > -10 ? prevBounds.z_fut : -8;
  }
  var prevHi = prevBounds.z_eff;

  var tCurr = infoFracs[k - 1]; /* info fraction at look k (0-indexed) */
  var tPrev = infoFracs[k - 2]; /* info fraction at look k-1 */
  var rho = Math.sqrt(tPrev / tCurr);
  var sigma = Math.sqrt(1 - tPrev / tCurr);

  return function(z) {
    /* f_k(z) = integral over z' in continuation region of
       f_{k-1}(z') * phi((z - rho*z') / sigma) / sigma dz' */
    return gaussLegendre(function(zp) {
      return prevDensity(zp) * normalPDF((z - rho * zp) / sigma) / sigma;
    }, prevLo, prevHi);
  };
}

/**
 * Compute P(Z_k > c | survived all previous looks) under H0.
 */
function computeCrossingProb(c, densityFn, prevLo, prevHi, rho, sigma, twoSided) {
  /* P = integral over z in [prevLo, prevHi] of
         densityFn(z) * P(Z_k > c | Z_{k-1} = z) dz
       = integral of densityFn(z) * (1 - Phi((c - rho*z)/sigma)) dz */
  var prob = gaussLegendre(function(z) {
    return densityFn(z) * (1 - normalCDF((c - rho * z) / sigma));
  }, prevLo, prevHi);

  if (twoSided) {
    /* Also count lower tail: P(Z_k < -c | survived) */
    prob += gaussLegendre(function(z) {
      return densityFn(z) * normalCDF((-c - rho * z) / sigma);
    }, prevLo, prevHi);
  }
  return prob;
}

/**
 * Compute P(Z_k < c | survived) under H0 for futility.
 */
function computeFutilityCrossingProb(c, densityFn, prevLo, prevHi, rho, sigma) {
  return gaussLegendre(function(z) {
    return densityFn(z) * normalCDF((c - rho * z) / sigma);
  }, prevLo, prevHi);
}


/* ===== 5. SVG CHART RENDERING ===== */


export { xoshiro128ss, seedPRNG, clamp, normalCDF, normalPDF, normalQuantile, normalRandom,
  alphaSpendOBF, alphaSpendPocock, alphaSpendLanDeMets, alphaSpendHSD, getSpendingFunction,
  gaussLegendre, computeBoundaries, buildRecursiveDensity, computeCrossingProb, computeFutilityCrossingProb };
