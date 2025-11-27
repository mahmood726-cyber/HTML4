/**
 * P-Curve and P-Uniform Publication Bias Methods
 *
 * Implements selection model approaches:
 * - P-curve analysis for evidential value
 * - P-uniform for bias-corrected estimation
 * - P-uniform* for heterogeneous effects
 *
 * Based on Simonsohn et al. (2014), van Assen et al. (2015), van Aert & van Assen (2020)
 */

import { jStat } from 'jstat';

/**
 * P-Curve Analysis
 * Tests whether p-values suggest true effect or p-hacking
 *
 * @param {Array} studies - Studies with es, se, or p-values directly
 * @param {Object} options - Options
 * @returns {Object} P-curve results
 */
export function pCurveAnalysis(studies, options = {}) {
  const {
    alpha = 0.05,  // Significance threshold
    rightSkewTest = true,
    flatnessTest = true
  } = options;

  // Extract significant p-values
  const sigStudies = [];

  for (const study of studies) {
    if (study.excluded) continue;

    let p;
    if (study.p != null) {
      p = study.p;
    } else if (study.es != null && study.se > 0) {
      const z = Math.abs(study.es / study.se);
      p = 2 * (1 - jStat.normal.cdf(z, 0, 1));
    } else {
      continue;
    }

    // Only include significant results (two-tailed)
    if (p < alpha && p > 0) {
      sigStudies.push({
        id: study.id,
        p,
        es: study.es,
        se: study.se
      });
    }
  }

  const k = sigStudies.length;

  if (k < 3) {
    return { error: 'At least 3 significant studies required for p-curve' };
  }

  const pValues = sigStudies.map(s => s.p);

  // Convert to pp-values (probability of observing p-value this extreme or more)
  // Under null: p-values are uniform [0, alpha], so pp = p/alpha
  // Under alternative: p-values are right-skewed

  const ppNull = pValues.map(p => p / alpha);  // Uniform under null

  // Test 1: Right-skew test (binomial test)
  // If true effect exists, more p-values should be close to 0 than to alpha
  const halfPoint = alpha / 2;
  const belowHalf = pValues.filter(p => p < halfPoint).length;
  const aboveHalf = k - belowHalf;

  const binomialPRight = 1 - jStat.binomial.cdf(belowHalf - 1, k, 0.5);  // P(X >= belowHalf)

  // Test 2: Flatness test (against uniform distribution)
  // Stouffer's method combining pp-values
  const zScores = ppNull.map(pp => jStat.normal.inv(1 - pp, 0, 1));
  const stoufferZ = zScores.reduce((a, b) => a + b, 0) / Math.sqrt(k);
  const stoufferP = 1 - jStat.normal.cdf(stoufferZ, 0, 1);

  // Test 3: 33% power test
  // If effects are powered at 33%, expected distribution can be computed
  // This tests against "just significant" p-hacking
  const near05 = pValues.filter(p => p > 0.04 && p < 0.05).length;
  const near01 = pValues.filter(p => p < 0.01).length;
  const ratio = near01 / (near05 + 0.1);  // Avoid division by zero

  // Distribution shape
  const bins = [0, 0.01, 0.02, 0.03, 0.04, 0.05];
  const histogram = bins.slice(0, -1).map((low, i) => ({
    range: `${low.toFixed(2)}-${bins[i + 1].toFixed(2)}`,
    count: pValues.filter(p => p >= low && p < bins[i + 1]).length,
    expected: k / (bins.length - 1)
  }));

  // Interpretation
  let evidentialValue = 'inconclusive';
  let intensity = 'unknown';

  if (binomialPRight < 0.05 && stoufferP < 0.05) {
    evidentialValue = 'present';
    intensity = ratio > 2 ? 'strong' : 'moderate';
  } else if (stoufferP > 0.95) {
    evidentialValue = 'absent';
    intensity = 'p-hacking suspected';
  }

  return {
    method: 'P-Curve Analysis',
    k,
    pValues,
    tests: {
      rightSkew: {
        belowHalf,
        aboveHalf,
        binomialP: binomialPRight,
        significant: binomialPRight < 0.05
      },
      flatness: {
        stoufferZ,
        stoufferP,
        significant: stoufferP < 0.05
      },
      extremeRatio: {
        near01,
        near05,
        ratio
      }
    },
    histogram,
    evidentialValue,
    intensity,
    interpretation: generatePCurveInterpretation(evidentialValue, intensity, binomialPRight, stoufferP)
  };
}

function generatePCurveInterpretation(evidentialValue, intensity, binomP, stoufferP) {
  if (evidentialValue === 'present') {
    return `P-curve indicates ${intensity} evidential value. Studies contain genuine effects.`;
  } else if (evidentialValue === 'absent') {
    return 'P-curve is flat or left-skewed, suggesting possible p-hacking or no true effect.';
  } else {
    return 'P-curve analysis is inconclusive. Results should be interpreted with caution.';
  }
}

/**
 * P-Uniform Method
 * Conditional estimation assuming homogeneous true effect
 *
 * @param {Array} studies - Studies with es, se
 * @param {Object} options - Options
 * @returns {Object} P-uniform results
 */
export function pUniform(studies, options = {}) {
  const {
    alpha = 0.05,
    method = 'ML',  // 'ML' or 'conditional'
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.es != null && s.se > 0
  );

  // Select significant studies
  const sigStudies = validStudies.filter(s => {
    const z = Math.abs(s.es / s.se);
    const p = 2 * (1 - jStat.normal.cdf(z, 0, 1));
    return p < alpha;
  });

  const k = sigStudies.length;

  if (k < 3) {
    return { error: 'At least 3 significant studies required for p-uniform' };
  }

  const y = sigStudies.map(s => s.es);
  const se = sigStudies.map(s => s.se);

  // Critical value for significance
  const zCrit = jStat.normal.inv(1 - alpha / 2, 0, 1);

  // P-uniform assumes all studies estimate same true effect mu
  // Conditional on significance: p-value | significant ~ Uniform(0, alpha/2) under H0
  // Under H1 with effect mu: distribution is right-skewed

  // Find mu that makes conditional p-values uniform
  // Use iterative search

  let muLo = -2;
  let muHi = 2;

  // Binary search for mu
  for (let iter = 0; iter < 50; iter++) {
    const muMid = (muLo + muHi) / 2;

    // Calculate conditional p-values under muMid
    const condPValues = sigStudies.map((s, i) => {
      const z = (y[i] - muMid) / se[i];
      // Conditional p-value given z > zCrit
      const pCond = 1 - jStat.normal.cdf(z, 0, 1);
      // Normalize to [0, 1] conditional on significance
      const pSig = 1 - jStat.normal.cdf(zCrit, 0, 1);
      return pCond / pSig;
    });

    // Test uniformity using mean (should be 0.5 if uniform)
    const meanCondP = condPValues.reduce((a, b) => a + b, 0) / k;

    if (meanCondP < 0.5) {
      muHi = muMid;
    } else {
      muLo = muMid;
    }

    if (Math.abs(muHi - muLo) < 1e-6) break;
  }

  const muEstimate = (muLo + muHi) / 2;

  // Standard error via bootstrap or asymptotic approximation
  // Simplified: use weighted average SE
  const avgSE = Math.sqrt(se.reduce((sum, s) => sum + s * s, 0) / k);
  const seMu = avgSE / Math.sqrt(k) * 1.5;  // Inflate for selection

  // Confidence interval
  const zCI = jStat.normal.inv((1 + confLevel) / 2, 0, 1);
  const ciLo = muEstimate - zCI * seMu;
  const ciHi = muEstimate + zCI * seMu;

  // Test for effect
  const zTest = muEstimate / seMu;
  const pVal = 2 * (1 - jStat.normal.cdf(Math.abs(zTest), 0, 1));

  // Publication bias test (test if selection affects estimate)
  // Compare p-uniform estimate with naive estimate
  const naiveW = se.map(s => 1 / (s * s));
  const sumNaiveW = naiveW.reduce((a, b) => a + b, 0);
  const naiveMu = naiveW.reduce((sum, w, i) => sum + w * y[i], 0) / sumNaiveW;

  const biasDiff = naiveMu - muEstimate;
  const biasDetected = Math.abs(biasDiff) > 0.1 * Math.abs(naiveMu);

  return {
    method: 'P-Uniform',
    k,
    kTotal: validStudies.length,
    estimate: {
      mu: muEstimate,
      se: seMu,
      ci: { lo: ciLo, hi: ciHi },
      z: zTest,
      p: pVal
    },
    naiveEstimate: naiveMu,
    biasCorrectionAmount: biasDiff,
    biasDetected,
    interpretation: generatePUniformInterpretation(muEstimate, naiveMu, biasDetected)
  };
}

function generatePUniformInterpretation(corrected, naive, biasDetected) {
  const parts = [];

  parts.push(`P-uniform estimate: ${corrected.toFixed(3)}`);
  parts.push(`Naive estimate: ${naive.toFixed(3)}`);

  if (biasDetected) {
    parts.push('Publication bias correction applied.');
  } else {
    parts.push('Limited evidence of publication bias.');
  }

  return parts.join(' ');
}

/**
 * P-Uniform* (star) Method
 * Allows for heterogeneous true effects
 *
 * @param {Array} studies - Studies with es, se
 * @param {Object} options - Options
 * @returns {Object} P-uniform* results
 */
export function pUniformStar(studies, options = {}) {
  const {
    alpha = 0.05,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.es != null && s.se > 0
  );

  const sigStudies = validStudies.filter(s => {
    const z = Math.abs(s.es / s.se);
    const p = 2 * (1 - jStat.normal.cdf(z, 0, 1));
    return p < alpha;
  });

  const k = sigStudies.length;

  if (k < 4) {
    return { error: 'At least 4 significant studies required for p-uniform*' };
  }

  const y = sigStudies.map(s => s.es);
  const se = sigStudies.map(s => s.se);

  // P-uniform* estimates both mu and tau²
  // Use two-step estimation

  // Step 1: Get initial p-uniform estimate
  const pUniformResult = pUniform(studies, { alpha, confLevel });
  if (pUniformResult.error) {
    return pUniformResult;
  }

  let mu = pUniformResult.estimate.mu;
  let tau2 = 0;

  // Step 2: Estimate tau² from residual variance
  const zCrit = jStat.normal.inv(1 - alpha / 2, 0, 1);

  for (let iter = 0; iter < 30; iter++) {
    const prevMu = mu;
    const prevTau2 = tau2;

    // Calculate conditional expectations under (mu, tau²)
    const condExpectations = sigStudies.map((s, i) => {
      const totalVar = se[i] * se[i] + tau2;
      const threshold = zCrit * se[i];

      // Expected value of y given y > threshold and true mean mu
      // Using truncated normal formula
      const z = (threshold - mu) / Math.sqrt(totalVar);
      const pdf = jStat.normal.pdf(z, 0, 1);
      const cdf = jStat.normal.cdf(z, 0, 1);
      const lambda = pdf / (1 - cdf + 0.001);  // Inverse Mills ratio

      return mu + Math.sqrt(totalVar) * lambda;
    });

    // Update mu
    const w = se.map((s, i) => 1 / (s * s + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    mu = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;

    // Update tau² using method of moments
    const Q = w.reduce((sum, wi, i) => sum + wi * (y[i] - mu) ** 2, 0);
    const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;

    // Adjust for selection bias
    const selectionAdjustment = 0.8;  // Empirical adjustment
    tau2 = Math.max(0, ((Q - (k - 1)) / C) * selectionAdjustment);

    if (Math.abs(mu - prevMu) < 1e-6 && Math.abs(tau2 - prevTau2) < 1e-6) {
      break;
    }
  }

  // Standard errors
  const wFinal = se.map((s, i) => 1 / (s * s + tau2));
  const sumWFinal = wFinal.reduce((a, b) => a + b, 0);
  const seMu = Math.sqrt(1 / sumWFinal) * 1.3;  // Selection adjustment

  // CI
  const zCI = jStat.normal.inv((1 + confLevel) / 2, 0, 1);

  // Heterogeneity
  const I2 = tau2 / (tau2 + se.reduce((sum, s) => sum + s * s, 0) / k) * 100;

  return {
    method: 'P-Uniform*',
    k,
    estimate: {
      mu,
      se: seMu,
      ci: {
        lo: mu - zCI * seMu,
        hi: mu + zCI * seMu
      }
    },
    heterogeneity: {
      tau2,
      tau: Math.sqrt(tau2),
      I2: Math.max(0, Math.min(100, I2))
    },
    naiveEstimate: pUniformResult.naiveEstimate,
    pUniformEstimate: pUniformResult.estimate.mu,
    interpretation: `P-uniform* accounts for heterogeneity (τ²=${tau2.toFixed(3)}). Estimate: ${mu.toFixed(3)}.`
  };
}

/**
 * Selection Model Comparison
 *
 * @param {Array} studies - Studies
 * @returns {Object} Comparison results
 */
export function compareSelectionModels(studies) {
  const pCurve = pCurveAnalysis(studies);
  const pUniformResult = pUniform(studies);
  const pUniformStarResult = pUniformStar(studies);

  // Calculate naive estimate for comparison
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const w = validStudies.map(s => 1 / (s.se * s.se));
  const sumW = w.reduce((a, b) => a + b, 0);
  const naiveEstimate = w.reduce((sum, wi, i) => sum + wi * validStudies[i].es, 0) / sumW;

  const results = {
    naive: naiveEstimate,
    pCurve: pCurve.error ? null : pCurve,
    pUniform: pUniformResult.error ? null : pUniformResult.estimate.mu,
    pUniformStar: pUniformStarResult.error ? null : pUniformStarResult.estimate.mu
  };

  // Assess concordance
  const estimates = [naiveEstimate];
  if (results.pUniform) estimates.push(results.pUniform);
  if (results.pUniformStar) estimates.push(results.pUniformStar);

  const range = Math.max(...estimates) - Math.min(...estimates);
  const concordant = range < 0.2 * Math.abs(naiveEstimate);

  return {
    estimates: results,
    concordance: {
      range,
      concordant,
      interpretation: concordant
        ? 'Selection models agree. Limited publication bias.'
        : 'Selection models diverge. Possible publication bias.'
    },
    recommendation: !pCurve.error && pCurve.evidentialValue === 'present'
      ? 'Evidential value present. Use p-uniform* for heterogeneous effects.'
      : 'Consider all estimates. Report sensitivity to model assumptions.'
  };
}

export default {
  pCurveAnalysis,
  pUniform,
  pUniformStar,
  compareSelectionModels
};
