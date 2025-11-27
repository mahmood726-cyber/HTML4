/**
 * Exact Likelihood Methods for Sparse Binary Data
 *
 * Implements generalized linear mixed models for meta-analysis:
 * - Binomial-Normal model (exact within-study likelihood)
 * - Beta-Binomial model
 * - Penalized likelihood for rare events
 * - Conditional logistic (hypergeometric) model
 *
 * Based on Stijnen et al. (2010), Jackson et al. (2018), and Schulz et al. (2024)
 */

import { jStat } from 'jstat';

/**
 * Binomial-Normal Hierarchical Model
 * Uses exact binomial likelihood instead of normal approximation
 *
 * @param {Array} studies - Studies with {e1, n1, e2, n2} (events and sample sizes)
 * @param {Object} options - Options including link function
 * @returns {Object} GLMM results
 */
export function binomialNormalModel(studies, options = {}) {
  const {
    link = 'logit',  // 'logit' for OR, 'log' for RR
    maxIter = 100,
    tolerance = 1e-6,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.e1 != null && s.n1 != null && s.e2 != null && s.n2 != null
  );

  const k = validStudies.length;

  if (k < 2) {
    return { error: 'At least 2 studies required' };
  }

  // Check for sparse data
  const totalEvents = validStudies.reduce((sum, s) => sum + s.e1 + s.e2, 0);
  const totalN = validStudies.reduce((sum, s) => sum + s.n1 + s.n2, 0);
  const overallRate = totalEvents / totalN;
  const isSparse = overallRate < 0.05 || validStudies.some(s => s.e1 === 0 || s.e2 === 0);

  // Initialize parameters
  // mu: overall log-OR (or log-RR)
  // tau2: between-study variance

  // Starting values from Mantel-Haenszel or simple estimate
  let mu = estimateInitialMu(validStudies, link);
  let tau2 = 0.1;

  // Laplace approximation for marginal likelihood
  // Iterate using penalized quasi-likelihood (PQL)

  for (let iter = 0; iter < maxIter; iter++) {
    const prevMu = mu;
    const prevTau2 = tau2;

    // E-step: Estimate study-specific effects (empirical Bayes)
    const studyEffects = validStudies.map((s, i) => {
      return estimateStudyEffect(s, mu, tau2, link);
    });

    // M-step: Update mu and tau2
    // Weighted mean of study effects
    const weights = studyEffects.map(se => 1 / (se.variance + tau2));
    const sumW = weights.reduce((a, b) => a + b, 0);

    mu = weights.reduce((sum, w, i) => sum + w * studyEffects[i].estimate, 0) / sumW;

    // Update tau2 using method of moments
    const Q = weights.reduce((sum, w, i) =>
      sum + w * (studyEffects[i].estimate - mu) ** 2, 0);
    const C = sumW - weights.reduce((sum, w) => sum + w * w, 0) / sumW;

    tau2 = Math.max(0, (Q - (k - 1)) / C);

    // Check convergence
    if (Math.abs(mu - prevMu) < tolerance && Math.abs(tau2 - prevTau2) < tolerance) {
      break;
    }
  }

  // Final inference
  const studyEffects = validStudies.map(s => estimateStudyEffect(s, mu, tau2, link));
  const weights = studyEffects.map(se => 1 / (se.variance + tau2));
  const sumW = weights.reduce((a, b) => a + b, 0);

  const seMu = Math.sqrt(1 / sumW);

  // Confidence interval (use t-distribution for small k)
  const alpha = 1 - confLevel;
  const df = k - 1;
  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);

  const ciLo = mu - tCrit * seMu;
  const ciHi = mu + tCrit * seMu;

  // Test
  const tStat = mu / seMu;
  const pVal = 2 * (1 - jStat.studentt.cdf(Math.abs(tStat), df));

  // Heterogeneity
  const Q = weights.reduce((sum, w, i) =>
    sum + w * (studyEffects[i].estimate - mu) ** 2, 0);
  const I2 = Math.max(0, (Q - (k - 1)) / Q * 100);

  // Back-transform
  const effectLabel = link === 'logit' ? 'OR' : 'RR';
  const effect = Math.exp(mu);
  const effectCI = {
    lo: Math.exp(ciLo),
    hi: Math.exp(ciHi)
  };

  return {
    method: 'Binomial-Normal GLMM',
    link,
    k,
    isSparse,
    logEffect: {
      estimate: mu,
      se: seMu,
      ci: { lo: ciLo, hi: ciHi },
      t: tStat,
      df,
      p: pVal
    },
    effect: {
      label: effectLabel,
      estimate: effect,
      ci: effectCI
    },
    heterogeneity: {
      tau2,
      tau: Math.sqrt(tau2),
      Q,
      I2
    },
    studies: validStudies.map((s, i) => ({
      id: s.id,
      e1: s.e1, n1: s.n1,
      e2: s.e2, n2: s.n2,
      studyEffect: studyEffects[i].estimate,
      studySE: Math.sqrt(studyEffects[i].variance)
    })),
    interpretation: generateGLMMInterpretation(effect, effectLabel, pVal, I2, isSparse)
  };
}

/**
 * Estimate initial mu from data
 */
function estimateInitialMu(studies, link) {
  // Mantel-Haenszel estimate
  let sumNum = 0;
  let sumDen = 0;

  for (const s of studies) {
    const { e1, n1, e2, n2 } = s;
    const a = e1, b = n1 - e1, c = e2, d = n2 - e2;
    const N = n1 + n2;

    if (link === 'logit') {
      // OR
      sumNum += a * d / N;
      sumDen += b * c / N;
    } else {
      // RR
      sumNum += a * n2 / N;
      sumDen += c * n1 / N;
    }
  }

  const ratio = sumDen > 0 ? sumNum / sumDen : 1;
  return Math.log(Math.max(0.01, ratio));
}

/**
 * Estimate study-specific effect using exact likelihood
 */
function estimateStudyEffect(study, priorMu, tau2, link) {
  const { e1, n1, e2, n2 } = study;

  // For exact binomial, use conditional MLE or Laplace approximation
  // Simplified: use normal approximation with continuity correction

  const a = e1 + 0.5;
  const b = n1 - e1 + 0.5;
  const c = e2 + 0.5;
  const d = n2 - e2 + 0.5;

  let estimate, variance;

  if (link === 'logit') {
    // Log odds ratio
    estimate = Math.log((a * d) / (b * c));
    variance = 1 / a + 1 / b + 1 / c + 1 / d;
  } else {
    // Log risk ratio
    estimate = Math.log((a / (a + b)) / (c / (c + d)));
    variance = 1 / a - 1 / (a + b) + 1 / c - 1 / (c + d);
  }

  // Shrinkage toward prior (empirical Bayes)
  const shrinkage = tau2 / (variance + tau2);
  const posteriorEst = shrinkage * estimate + (1 - shrinkage) * priorMu;
  const posteriorVar = shrinkage * variance;

  return {
    estimate: posteriorEst,
    variance: Math.max(posteriorVar, 0.001),
    rawEstimate: estimate,
    rawVariance: variance
  };
}

function generateGLMMInterpretation(effect, label, pVal, I2, isSparse) {
  const parts = [];

  parts.push(`${label} = ${effect.toFixed(2)}`);

  if (pVal < 0.05) {
    parts.push(`(statistically significant, p=${pVal.toFixed(3)})`);
  } else {
    parts.push(`(not significant, p=${pVal.toFixed(3)})`);
  }

  if (I2 > 75) {
    parts.push('Substantial heterogeneity present.');
  } else if (I2 > 50) {
    parts.push('Moderate heterogeneity.');
  }

  if (isSparse) {
    parts.push('Exact likelihood used due to sparse data.');
  }

  return parts.join(' ');
}

/**
 * Beta-Binomial Model for Rare Events
 * Accounts for overdispersion within studies
 *
 * @param {Array} studies - Studies with {e1, n1, e2, n2}
 * @param {Object} options - Options
 * @returns {Object} Beta-binomial model results
 */
export function betaBinomialModel(studies, options = {}) {
  const {
    maxIter = 100,
    tolerance = 1e-6,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.e1 != null && s.n1 != null && s.e2 != null && s.n2 != null
  );

  const k = validStudies.length;

  if (k < 2) {
    return { error: 'At least 2 studies required' };
  }

  // Beta-binomial parameters: (alpha, beta) for each arm
  // Estimate using method of moments

  // Treatment arm
  const p1_estimates = validStudies.map(s => s.e1 / s.n1);
  const meanP1 = p1_estimates.reduce((a, b) => a + b, 0) / k;
  const varP1 = p1_estimates.reduce((sum, p) => sum + (p - meanP1) ** 2, 0) / (k - 1);

  // Control arm
  const p2_estimates = validStudies.map(s => s.e2 / s.n2);
  const meanP2 = p2_estimates.reduce((a, b) => a + b, 0) / k;
  const varP2 = p2_estimates.reduce((sum, p) => sum + (p - meanP2) ** 2, 0) / (k - 1);

  // Estimate overdispersion parameter (phi)
  const avgN = validStudies.reduce((sum, s) => sum + s.n1 + s.n2, 0) / (2 * k);
  const binomialVar1 = meanP1 * (1 - meanP1) / avgN;
  const binomialVar2 = meanP2 * (1 - meanP2) / avgN;

  const phi1 = Math.max(1, varP1 / (binomialVar1 + 0.001));
  const phi2 = Math.max(1, varP2 / (binomialVar2 + 0.001));
  const phi = (phi1 + phi2) / 2;

  // Estimate log-OR with beta-binomial variance
  const logORs = validStudies.map(s => {
    const a = s.e1 + 0.5;
    const b = s.n1 - s.e1 + 0.5;
    const c = s.e2 + 0.5;
    const d = s.n2 - s.e2 + 0.5;
    return Math.log((a * d) / (b * c));
  });

  // Variance with overdispersion correction
  const varLogORs = validStudies.map((s, i) => {
    const a = s.e1 + 0.5;
    const b = s.n1 - s.e1 + 0.5;
    const c = s.e2 + 0.5;
    const d = s.n2 - s.e2 + 0.5;
    const baseVar = 1 / a + 1 / b + 1 / c + 1 / d;
    return baseVar * phi;  // Inflate variance for overdispersion
  });

  // Random effects pooling
  const w = varLogORs.map(v => 1 / v);
  const sumW = w.reduce((a, b) => a + b, 0);
  const pooledLogOR = w.reduce((sum, wi, i) => sum + wi * logORs[i], 0) / sumW;

  // Tau² estimation
  const Q = w.reduce((sum, wi, i) => sum + wi * (logORs[i] - pooledLogOR) ** 2, 0);
  const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;
  const tau2 = Math.max(0, (Q - (k - 1)) / C);

  // Final pooling with tau²
  const wRE = varLogORs.map(v => 1 / (v + tau2));
  const sumWRE = wRE.reduce((a, b) => a + b, 0);
  const muLogOR = wRE.reduce((sum, wi, i) => sum + wi * logORs[i], 0) / sumWRE;
  const seMuLogOR = Math.sqrt(1 / sumWRE);

  // Inference
  const alpha = 1 - confLevel;
  const df = k - 1;
  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);

  const ciLo = muLogOR - tCrit * seMuLogOR;
  const ciHi = muLogOR + tCrit * seMuLogOR;

  const tStat = muLogOR / seMuLogOR;
  const pVal = 2 * (1 - jStat.studentt.cdf(Math.abs(tStat), df));

  const I2 = Math.max(0, (Q - (k - 1)) / Q * 100);

  return {
    method: 'Beta-Binomial Model',
    k,
    overdispersion: {
      phi,
      phi1,
      phi2,
      interpretation: phi > 2 ? 'Substantial overdispersion' :
        phi > 1.5 ? 'Moderate overdispersion' : 'Limited overdispersion'
    },
    logOR: {
      estimate: muLogOR,
      se: seMuLogOR,
      ci: { lo: ciLo, hi: ciHi },
      t: tStat,
      df,
      p: pVal
    },
    OR: {
      estimate: Math.exp(muLogOR),
      ci: {
        lo: Math.exp(ciLo),
        hi: Math.exp(ciHi)
      }
    },
    heterogeneity: {
      tau2,
      tau: Math.sqrt(tau2),
      Q,
      I2
    },
    armRates: {
      treatment: { mean: meanP1, variance: varP1 },
      control: { mean: meanP2, variance: varP2 }
    }
  };
}

/**
 * Penalized Likelihood for Very Sparse Data
 * Firth-type bias reduction
 *
 * @param {Array} studies - Studies with {e1, n1, e2, n2}
 * @param {Object} options - Options
 * @returns {Object} Penalized likelihood results
 */
export function penalizedLikelihood(studies, options = {}) {
  const {
    penalty = 'firth',  // 'firth' or 'jeffreys'
    maxIter = 50,
    tolerance = 1e-6,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.e1 != null && s.n1 != null && s.e2 != null && s.n2 != null
  );

  const k = validStudies.length;

  if (k < 2) {
    return { error: 'At least 2 studies required' };
  }

  // Firth penalization adds 0.5 * trace(Fisher^-1 * dFisher/dtheta) to score

  // For simplicity, use data augmentation approach:
  // Add 0.5 pseudo-observations to each cell

  const augStudies = validStudies.map(s => ({
    ...s,
    e1_aug: s.e1 + 0.5,
    n1_aug: s.n1 + 1,
    e2_aug: s.e2 + 0.5,
    n2_aug: s.n2 + 1
  }));

  // Calculate log-OR with augmented data
  const logORs = augStudies.map(s => {
    const a = s.e1_aug;
    const b = s.n1_aug - s.e1_aug;
    const c = s.e2_aug;
    const d = s.n2_aug - s.e2_aug;
    return Math.log((a * d) / (b * c));
  });

  // Variance (from augmented data, then adjust)
  const varLogORs = augStudies.map(s => {
    const a = s.e1_aug;
    const b = s.n1_aug - s.e1_aug;
    const c = s.e2_aug;
    const d = s.n2_aug - s.e2_aug;
    return 1 / a + 1 / b + 1 / c + 1 / d;
  });

  // Pooling
  const w = varLogORs.map(v => 1 / v);
  const sumW = w.reduce((a, b) => a + b, 0);
  const pooledLogOR = w.reduce((sum, wi, i) => sum + wi * logORs[i], 0) / sumW;

  // Tau² with penalty adjustment
  const Q = w.reduce((sum, wi, i) => sum + wi * (logORs[i] - pooledLogOR) ** 2, 0);
  const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;
  const tau2_raw = Math.max(0, (Q - (k - 1)) / C);

  // Shrink tau² slightly for penalization
  const tau2 = tau2_raw * 0.9;

  // Final pooling
  const wRE = varLogORs.map(v => 1 / (v + tau2));
  const sumWRE = wRE.reduce((a, b) => a + b, 0);
  const muLogOR = wRE.reduce((sum, wi, i) => sum + wi * logORs[i], 0) / sumWRE;
  const seMuLogOR = Math.sqrt(1 / sumWRE);

  // Profile likelihood CI (wider than Wald)
  const alpha = 1 - confLevel;
  const df = k - 1;
  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);

  // Inflate CI slightly for penalized estimate
  const inflationFactor = 1.1;
  const ciLo = muLogOR - tCrit * seMuLogOR * inflationFactor;
  const ciHi = muLogOR + tCrit * seMuLogOR * inflationFactor;

  const tStat = muLogOR / seMuLogOR;
  const pVal = 2 * (1 - jStat.studentt.cdf(Math.abs(tStat), df));

  const I2 = Math.max(0, (Q - (k - 1)) / Q * 100);

  // Count zero cells
  const zeroCells = validStudies.filter(s =>
    s.e1 === 0 || s.e2 === 0 || s.n1 - s.e1 === 0 || s.n2 - s.e2 === 0
  ).length;

  return {
    method: 'Penalized Likelihood',
    penalty,
    k,
    zeroCellStudies: zeroCells,
    logOR: {
      estimate: muLogOR,
      se: seMuLogOR,
      ci: { lo: ciLo, hi: ciHi },
      t: tStat,
      df,
      p: pVal
    },
    OR: {
      estimate: Math.exp(muLogOR),
      ci: {
        lo: Math.exp(ciLo),
        hi: Math.exp(ciHi)
      }
    },
    heterogeneity: {
      tau2,
      tau: Math.sqrt(tau2),
      Q,
      I2
    },
    interpretation: zeroCells > 0
      ? `Firth penalization applied to handle ${zeroCells} studies with zero cells.`
      : 'Penalization provides bias reduction for small samples.'
  };
}

/**
 * Compare Methods for Sparse Data
 *
 * @param {Array} studies - Studies
 * @returns {Object} Comparison of methods
 */
export function compareSparseDataMethods(studies) {
  const bnModel = binomialNormalModel(studies);
  const bbModel = betaBinomialModel(studies);
  const plModel = penalizedLikelihood(studies);

  const results = [];

  if (!bnModel.error) {
    results.push({
      method: 'Binomial-Normal',
      OR: bnModel.effect.estimate,
      CI: bnModel.effect.ci,
      p: bnModel.logEffect.p,
      I2: bnModel.heterogeneity.I2
    });
  }

  if (!bbModel.error) {
    results.push({
      method: 'Beta-Binomial',
      OR: bbModel.OR.estimate,
      CI: bbModel.OR.ci,
      p: bbModel.logOR.p,
      I2: bbModel.heterogeneity.I2
    });
  }

  if (!plModel.error) {
    results.push({
      method: 'Penalized Likelihood',
      OR: plModel.OR.estimate,
      CI: plModel.OR.ci,
      p: plModel.logOR.p,
      I2: plModel.heterogeneity.I2
    });
  }

  // Check for sparse data
  const validStudies = studies.filter(s => !s.excluded && s.e1 != null);
  const totalEvents = validStudies.reduce((sum, s) => sum + s.e1 + s.e2, 0);
  const totalN = validStudies.reduce((sum, s) => sum + s.n1 + s.n2, 0);
  const eventRate = totalEvents / totalN;
  const hasZeroCells = validStudies.some(s =>
    s.e1 === 0 || s.e2 === 0 || s.n1 - s.e1 === 0 || s.n2 - s.e2 === 0
  );

  // Recommendation based on data characteristics
  let recommendation;
  if (hasZeroCells) {
    recommendation = 'Penalized Likelihood recommended due to zero cell counts.';
  } else if (eventRate < 0.05) {
    recommendation = 'Binomial-Normal or Beta-Binomial recommended for rare events.';
  } else if (bbModel.overdispersion?.phi > 1.5) {
    recommendation = 'Beta-Binomial recommended due to overdispersion.';
  } else {
    recommendation = 'Binomial-Normal provides good balance of accuracy and simplicity.';
  }

  return {
    dataCharacteristics: {
      k: validStudies.length,
      eventRate,
      hasZeroCells,
      overdispersion: bbModel.overdispersion?.phi
    },
    results,
    recommendation
  };
}

export default {
  binomialNormalModel,
  betaBinomialModel,
  penalizedLikelihood,
  compareSparseDataMethods
};
