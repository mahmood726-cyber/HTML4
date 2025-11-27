/**
 * Dose-Response Meta-Analysis
 *
 * Methods for analyzing dose-response relationships:
 * - Greenland & Longnecker method for trend estimation
 * - Restricted cubic splines for non-linear dose-response
 * - One-stage and two-stage approaches
 * - Model comparison (linear, quadratic, spline)
 *
 * Based on dosresmeta methodology and GRADE guidance for dose-response
 */

import { jStat } from 'jstat';

/**
 * Calculate covariance matrix for dose categories within a study
 * Using Greenland & Longnecker (1992) method
 *
 * @param {Array} categories - Dose categories [{dose, cases, n, rr, ciLo, ciHi}]
 * @param {number} referenceIdx - Index of reference category (usually 0)
 * @returns {Object} Covariance matrix and log RRs
 */
export function calculateDoseCovariance(categories, referenceIdx = 0) {
  const k = categories.length;

  if (k < 2) {
    return { error: 'At least 2 dose categories required' };
  }

  // Extract data
  const doses = categories.map(c => c.dose);
  const cases = categories.map(c => c.cases);
  const n = categories.map(c => c.n);

  // Log RRs relative to reference
  const logRR = categories.map((c, i) => {
    if (i === referenceIdx) return 0;
    return c.rr ? Math.log(c.rr) : 0;
  });

  // Variances from CIs
  const varLogRR = categories.map((c, i) => {
    if (i === referenceIdx) return 0;
    if (c.ciLo && c.ciHi) {
      const seLogRR = (Math.log(c.ciHi) - Math.log(c.ciLo)) / (2 * 1.96);
      return seLogRR * seLogRR;
    }
    // Estimate from cases/n if no CI
    if (c.cases && c.n && c.cases > 0) {
      return 1 / c.cases + 1 / (c.n - c.cases);
    }
    return 0.1;  // Default
  });

  // Covariance matrix (Greenland & Longnecker)
  // Cov(log RRi, log RRj) = 1/n_ref for i ≠ j
  const refCases = cases[referenceIdx];
  const refN = n[referenceIdx];
  const refNonCases = refN - refCases;

  const covMatrix = [];
  for (let i = 0; i < k; i++) {
    covMatrix[i] = [];
    for (let j = 0; j < k; j++) {
      if (i === referenceIdx || j === referenceIdx) {
        covMatrix[i][j] = 0;
      } else if (i === j) {
        covMatrix[i][j] = varLogRR[i];
      } else {
        // Covariance due to shared reference
        covMatrix[i][j] = refCases > 0 ? 1 / refCases : 0;
      }
    }
  }

  return {
    doses,
    logRR,
    varLogRR,
    covMatrix,
    referenceIdx,
    referenceDose: doses[referenceIdx]
  };
}

/**
 * Linear dose-response meta-regression (one-stage)
 * Estimates trend per unit dose increase
 *
 * @param {Array} studies - Studies with dose categories
 * @param {Object} options - Options
 * @returns {Object} Linear trend estimate
 */
export function linearDoseResponse(studies, options = {}) {
  const { doseUnit = 1, referenceIdx = 0 } = options;

  const validStudies = studies.filter(s => s.categories && s.categories.length >= 2);

  if (validStudies.length < 2) {
    return { error: 'At least 2 studies with dose categories required' };
  }

  // Process each study
  const studyResults = [];

  for (const study of validStudies) {
    const covData = calculateDoseCovariance(study.categories, referenceIdx);
    if (covData.error) continue;

    // Exclude reference category
    const nonRefIdx = covData.doses.map((_, i) => i).filter(i => i !== referenceIdx);

    if (nonRefIdx.length === 0) continue;

    // Center doses
    const refDose = covData.referenceDose;
    const centeredDoses = nonRefIdx.map(i => (covData.doses[i] - refDose) / doseUnit);
    const y = nonRefIdx.map(i => covData.logRR[i]);
    const v = nonRefIdx.map(i => covData.varLogRR[i]);

    // Weighted regression through origin (log RR = β * dose)
    // β = Σ(wi * xi * yi) / Σ(wi * xi²)
    const w = v.map(vi => vi > 0 ? 1 / vi : 0);

    let sumWXY = 0, sumWX2 = 0;
    for (let i = 0; i < centeredDoses.length; i++) {
      sumWXY += w[i] * centeredDoses[i] * y[i];
      sumWX2 += w[i] * centeredDoses[i] * centeredDoses[i];
    }

    const beta = sumWX2 > 0 ? sumWXY / sumWX2 : 0;
    const seBeta = sumWX2 > 0 ? Math.sqrt(1 / sumWX2) : 0;

    studyResults.push({
      id: study.id,
      beta,
      seBeta,
      weight: sumWX2
    });
  }

  if (studyResults.length === 0) {
    return { error: 'No valid studies for dose-response analysis' };
  }

  // Pool study-level slopes (random effects)
  const betas = studyResults.map(s => s.beta);
  const vars = studyResults.map(s => s.seBeta * s.seBeta);
  const weights = vars.map(v => v > 0 ? 1 / v : 0);

  const sumW = weights.reduce((a, b) => a + b, 0);
  const pooledBeta = weights.reduce((sum, w, i) => sum + w * betas[i], 0) / sumW;

  // Heterogeneity
  const Q = weights.reduce((sum, w, i) => sum + w * (betas[i] - pooledBeta) ** 2, 0);
  const df = studyResults.length - 1;
  const pHet = 1 - jStat.chisquare.cdf(Q, df);

  // DerSimonian-Laird tau²
  const C = sumW - weights.reduce((sum, w) => sum + w * w, 0) / sumW;
  const tau2 = Math.max(0, (Q - df) / C);

  // Random effects pooled estimate
  const reWeights = vars.map(v => 1 / (v + tau2));
  const sumREW = reWeights.reduce((a, b) => a + b, 0);
  const rePooledBeta = reWeights.reduce((sum, w, i) => sum + w * betas[i], 0) / sumREW;
  const sePooledBeta = Math.sqrt(1 / sumREW);

  const I2 = df > 0 ? Math.max(0, (Q - df) / Q * 100) : 0;

  return {
    method: 'Linear Dose-Response',
    k: studyResults.length,
    doseUnit,
    fixedEffect: {
      beta: pooledBeta,
      se: Math.sqrt(1 / sumW),
      ci: {
        lo: pooledBeta - 1.96 * Math.sqrt(1 / sumW),
        hi: pooledBeta + 1.96 * Math.sqrt(1 / sumW)
      },
      rr: Math.exp(pooledBeta),
      rrCI: {
        lo: Math.exp(pooledBeta - 1.96 * Math.sqrt(1 / sumW)),
        hi: Math.exp(pooledBeta + 1.96 * Math.sqrt(1 / sumW))
      }
    },
    randomEffects: {
      beta: rePooledBeta,
      se: sePooledBeta,
      ci: {
        lo: rePooledBeta - 1.96 * sePooledBeta,
        hi: rePooledBeta + 1.96 * sePooledBeta
      },
      rr: Math.exp(rePooledBeta),
      rrCI: {
        lo: Math.exp(rePooledBeta - 1.96 * sePooledBeta),
        hi: Math.exp(rePooledBeta + 1.96 * sePooledBeta)
      }
    },
    heterogeneity: {
      Q,
      df,
      pHet,
      tau2,
      I2
    },
    studies: studyResults,
    interpretation: generateLinearInterpretation(rePooledBeta, sePooledBeta, doseUnit)
  };
}

function generateLinearInterpretation(beta, se, doseUnit) {
  const rr = Math.exp(beta);
  const significant = Math.abs(beta) > 1.96 * se;

  if (!significant) {
    return `No significant dose-response trend detected (RR per ${doseUnit} unit = ${rr.toFixed(2)})`;
  }

  if (beta > 0) {
    return `Positive dose-response: RR increases by ${((rr - 1) * 100).toFixed(1)}% per ${doseUnit} unit increase`;
  } else {
    return `Inverse dose-response: RR decreases by ${((1 - rr) * 100).toFixed(1)}% per ${doseUnit} unit increase`;
  }
}

/**
 * Quadratic dose-response model
 *
 * @param {Array} studies - Studies with dose categories
 * @param {Object} options - Options
 * @returns {Object} Quadratic model results
 */
export function quadraticDoseResponse(studies, options = {}) {
  const { doseUnit = 1, referenceIdx = 0 } = options;

  const validStudies = studies.filter(s => s.categories && s.categories.length >= 3);

  if (validStudies.length < 2) {
    return { error: 'At least 2 studies with 3+ dose categories required for quadratic model' };
  }

  // Collect all data points
  let allDoses = [];
  let allLogRR = [];
  let allWeights = [];
  let allStudyIdx = [];

  validStudies.forEach((study, sIdx) => {
    const covData = calculateDoseCovariance(study.categories, referenceIdx);
    if (covData.error) return;

    const refDose = covData.referenceDose;

    covData.doses.forEach((dose, i) => {
      if (i === referenceIdx) return;

      const centeredDose = (dose - refDose) / doseUnit;
      const w = covData.varLogRR[i] > 0 ? 1 / covData.varLogRR[i] : 0;

      allDoses.push(centeredDose);
      allLogRR.push(covData.logRR[i]);
      allWeights.push(w);
      allStudyIdx.push(sIdx);
    });
  });

  const n = allDoses.length;
  if (n < 4) {
    return { error: 'Insufficient data points for quadratic model' };
  }

  // Weighted least squares: y = β1*x + β2*x²
  // Normal equations
  let sumW = 0, sumWX = 0, sumWX2 = 0, sumWX3 = 0, sumWX4 = 0;
  let sumWY = 0, sumWXY = 0, sumWX2Y = 0;

  for (let i = 0; i < n; i++) {
    const w = allWeights[i];
    const x = allDoses[i];
    const x2 = x * x;
    const y = allLogRR[i];

    sumW += w;
    sumWX += w * x;
    sumWX2 += w * x2;
    sumWX3 += w * x * x2;
    sumWX4 += w * x2 * x2;
    sumWY += w * y;
    sumWXY += w * x * y;
    sumWX2Y += w * x2 * y;
  }

  // Solve 2x2 system for β1 and β2
  const det = sumWX2 * sumWX4 - sumWX3 * sumWX3;

  if (Math.abs(det) < 1e-10) {
    return { error: 'Singular matrix in quadratic regression' };
  }

  const beta1 = (sumWX4 * sumWXY - sumWX3 * sumWX2Y) / det;
  const beta2 = (sumWX2 * sumWX2Y - sumWX3 * sumWXY) / det;

  // Residual variance and SEs (simplified)
  let SSres = 0;
  for (let i = 0; i < n; i++) {
    const fitted = beta1 * allDoses[i] + beta2 * allDoses[i] * allDoses[i];
    SSres += allWeights[i] * (allLogRR[i] - fitted) ** 2;
  }

  const MSres = SSres / (n - 2);
  const seBeta1 = Math.sqrt(MSres * sumWX4 / det);
  const seBeta2 = Math.sqrt(MSres * sumWX2 / det);

  // Test for non-linearity (is β2 significant?)
  const tBeta2 = beta2 / seBeta2;
  const pNonLinear = 2 * (1 - jStat.studentt.cdf(Math.abs(tBeta2), n - 2));

  // Find minimum/maximum dose (for concave/convex curves)
  let extremum = null;
  if (Math.abs(beta2) > 0.001) {
    extremum = -beta1 / (2 * beta2);
  }

  // Generate curve for plotting
  const doseRange = [Math.min(...allDoses), Math.max(...allDoses)];
  const curve = [];
  for (let d = doseRange[0]; d <= doseRange[1]; d += (doseRange[1] - doseRange[0]) / 50) {
    const logRR = beta1 * d + beta2 * d * d;
    curve.push({
      dose: d * doseUnit,
      logRR,
      rr: Math.exp(logRR)
    });
  }

  return {
    method: 'Quadratic Dose-Response',
    nPoints: n,
    nStudies: validStudies.length,
    doseUnit,
    coefficients: {
      linear: {
        estimate: beta1,
        se: seBeta1,
        p: 2 * (1 - jStat.studentt.cdf(Math.abs(beta1 / seBeta1), n - 2))
      },
      quadratic: {
        estimate: beta2,
        se: seBeta2,
        p: pNonLinear
      }
    },
    nonLinearity: {
      pValue: pNonLinear,
      significant: pNonLinear < 0.05,
      shape: beta2 > 0 ? 'U-shaped' : 'inverted-U'
    },
    extremum: extremum ? {
      dose: extremum * doseUnit,
      type: beta2 > 0 ? 'minimum' : 'maximum',
      logRR: beta1 * extremum + beta2 * extremum * extremum,
      rr: Math.exp(beta1 * extremum + beta2 * extremum * extremum)
    } : null,
    curve,
    interpretation: generateQuadraticInterpretation(beta1, beta2, pNonLinear, doseUnit)
  };
}

function generateQuadraticInterpretation(beta1, beta2, pNonLinear, doseUnit) {
  if (pNonLinear >= 0.05) {
    return 'No significant non-linearity detected. Linear model may be adequate.';
  }

  if (beta2 > 0) {
    return `Significant U-shaped dose-response (p=${pNonLinear.toFixed(3)}). Risk decreases then increases with dose.`;
  } else {
    return `Significant inverted-U relationship (p=${pNonLinear.toFixed(3)}). Risk increases then decreases with dose.`;
  }
}

/**
 * Restricted Cubic Spline Dose-Response
 * For flexible non-linear modeling
 *
 * @param {Array} studies - Studies with dose categories
 * @param {Object} options - Options including knot positions
 * @returns {Object} Spline model results
 */
export function splineDoseResponse(studies, options = {}) {
  const {
    doseUnit = 1,
    referenceIdx = 0,
    nKnots = 3,  // Default 3 knots (percentiles 10, 50, 90)
    knotPositions = null
  } = options;

  // Collect data points
  let allDoses = [];
  let allLogRR = [];
  let allWeights = [];

  const validStudies = studies.filter(s => s.categories && s.categories.length >= 2);

  validStudies.forEach(study => {
    const covData = calculateDoseCovariance(study.categories, referenceIdx);
    if (covData.error) return;

    const refDose = covData.referenceDose;

    covData.doses.forEach((dose, i) => {
      if (i === referenceIdx) return;

      allDoses.push((dose - refDose) / doseUnit);
      allLogRR.push(covData.logRR[i]);
      allWeights.push(covData.varLogRR[i] > 0 ? 1 / covData.varLogRR[i] : 0);
    });
  });

  const n = allDoses.length;
  if (n < nKnots + 2) {
    return { error: `At least ${nKnots + 2} data points required for ${nKnots}-knot spline` };
  }

  // Determine knot positions
  const sortedDoses = [...allDoses].sort((a, b) => a - b);
  const knots = knotPositions || [
    jStat.percentile(sortedDoses, 10),
    jStat.percentile(sortedDoses, 50),
    jStat.percentile(sortedDoses, 90)
  ].slice(0, nKnots);

  // Create spline basis functions
  // RCS with 3 knots: 2 basis functions (x, and one truncated power function)
  function rcsTransform(x, knots) {
    const k = knots.length;
    const basis = [x];  // First basis is linear

    // Additional basis functions
    for (let j = 0; j < k - 2; j++) {
      const kj = knots[j];
      const kLast = knots[k - 1];
      const kPrev = knots[k - 2];

      const h = (Math.pow(Math.max(0, x - kj), 3) -
        Math.pow(Math.max(0, x - kPrev), 3) * (kLast - kj) / (kLast - kPrev) +
        Math.pow(Math.max(0, x - kLast), 3) * (kPrev - kj) / (kLast - kPrev)) /
        Math.pow(kLast - knots[0], 2);

      basis.push(h);
    }

    return basis;
  }

  // Build design matrix
  const X = allDoses.map(d => rcsTransform(d, knots));
  const nBasis = X[0].length;

  // Weighted least squares
  // (X'WX)β = X'Wy
  const XtWX = [];
  const XtWy = [];

  for (let i = 0; i < nBasis; i++) {
    XtWX[i] = [];
    XtWy[i] = 0;
    for (let j = 0; j < nBasis; j++) {
      XtWX[i][j] = 0;
      for (let k = 0; k < n; k++) {
        XtWX[i][j] += allWeights[k] * X[k][i] * X[k][j];
      }
    }
    for (let k = 0; k < n; k++) {
      XtWy[i] += allWeights[k] * X[k][i] * allLogRR[k];
    }
  }

  // Solve (simple 2x2 for 3 knots)
  let beta;
  if (nBasis === 2) {
    const det = XtWX[0][0] * XtWX[1][1] - XtWX[0][1] * XtWX[1][0];
    if (Math.abs(det) < 1e-10) {
      return { error: 'Singular matrix in spline regression' };
    }
    beta = [
      (XtWX[1][1] * XtWy[0] - XtWX[0][1] * XtWy[1]) / det,
      (XtWX[0][0] * XtWy[1] - XtWX[1][0] * XtWy[0]) / det
    ];
  } else {
    // For more complex cases, use simple iterative solution
    beta = new Array(nBasis).fill(0);
    // Gradient descent (simplified)
    for (let iter = 0; iter < 100; iter++) {
      const grad = XtWy.map((xwy, i) => {
        let sum = 0;
        for (let j = 0; j < nBasis; j++) {
          sum += XtWX[i][j] * beta[j];
        }
        return xwy - sum;
      });
      for (let i = 0; i < nBasis; i++) {
        beta[i] += 0.1 * grad[i] / (XtWX[i][i] + 0.001);
      }
    }
  }

  // Generate curve
  const doseRange = [Math.min(...allDoses), Math.max(...allDoses)];
  const curve = [];
  for (let d = doseRange[0]; d <= doseRange[1]; d += (doseRange[1] - doseRange[0]) / 100) {
    const basis = rcsTransform(d, knots);
    const logRR = basis.reduce((sum, b, i) => sum + b * beta[i], 0);
    curve.push({
      dose: d * doseUnit,
      logRR,
      rr: Math.exp(logRR)
    });
  }

  // Test for non-linearity (overall test)
  let SSlin = 0, SSspline = 0;
  for (let i = 0; i < n; i++) {
    const linFitted = beta[0] * allDoses[i];
    const splineFitted = X[i].reduce((sum, b, j) => sum + b * beta[j], 0);
    SSlin += allWeights[i] * (allLogRR[i] - linFitted) ** 2;
    SSspline += allWeights[i] * (allLogRR[i] - splineFitted) ** 2;
  }

  const FnonLin = ((SSlin - SSspline) / (nBasis - 1)) / (SSspline / (n - nBasis));
  const pNonLinear = 1 - jStat.centralF.cdf(FnonLin, nBasis - 1, n - nBasis);

  return {
    method: 'Restricted Cubic Spline',
    nPoints: n,
    nStudies: validStudies.length,
    doseUnit,
    nKnots,
    knots: knots.map(k => k * doseUnit),
    coefficients: beta.map((b, i) => ({
      basis: i,
      estimate: b
    })),
    nonLinearity: {
      Fstat: FnonLin,
      pValue: pNonLinear,
      significant: pNonLinear < 0.05
    },
    curve,
    interpretation: pNonLinear < 0.05
      ? 'Significant non-linear dose-response detected. Spline model preferred.'
      : 'No significant non-linearity. Linear model may be adequate.'
  };
}

/**
 * Compare dose-response models
 *
 * @param {Array} studies - Studies with dose categories
 * @param {Object} options - Options
 * @returns {Object} Model comparison results
 */
export function compareDoseResponseModels(studies, options = {}) {
  const linear = linearDoseResponse(studies, options);
  const quadratic = quadraticDoseResponse(studies, options);
  const spline = splineDoseResponse(studies, options);

  const results = [];

  if (!linear.error) {
    results.push({
      model: 'Linear',
      nParams: 1,
      interpretation: linear.interpretation
    });
  }

  if (!quadratic.error) {
    results.push({
      model: 'Quadratic',
      nParams: 2,
      pNonLinear: quadratic.nonLinearity.pValue,
      interpretation: quadratic.interpretation
    });
  }

  if (!spline.error) {
    results.push({
      model: 'Spline',
      nParams: spline.coefficients.length,
      pNonLinear: spline.nonLinearity.pValue,
      interpretation: spline.interpretation
    });
  }

  // Recommendation
  let recommended = 'Linear';
  if (!quadratic.error && quadratic.nonLinearity.significant) {
    recommended = 'Quadratic';
  }
  if (!spline.error && spline.nonLinearity.significant &&
    (!quadratic.error && spline.nonLinearity.pValue < quadratic.nonLinearity.pValue)) {
    recommended = 'Spline';
  }

  return {
    models: results,
    recommended,
    linear: linear.error ? { error: linear.error } : {
      beta: linear.randomEffects.beta,
      rr: linear.randomEffects.rr
    },
    quadratic: quadratic.error ? { error: quadratic.error } : {
      linear: quadratic.coefficients.linear.estimate,
      quadratic: quadratic.coefficients.quadratic.estimate,
      pNonLinear: quadratic.nonLinearity.pValue
    },
    spline: spline.error ? { error: spline.error } : {
      pNonLinear: spline.nonLinearity.pValue
    }
  };
}

/**
 * Generate dose-response plot data
 *
 * @param {Object} modelResult - Result from any dose-response model
 * @param {Array} studies - Original studies for data points
 * @returns {Object} Plot data
 */
export function getDoseResponsePlotData(modelResult, studies, options = {}) {
  const { referenceIdx = 0 } = options;

  // Study-level points
  const points = [];

  studies.forEach(study => {
    if (!study.categories) return;

    const covData = calculateDoseCovariance(study.categories, referenceIdx);
    if (covData.error) return;

    covData.doses.forEach((dose, i) => {
      if (i === referenceIdx) return;

      points.push({
        studyId: study.id,
        dose,
        logRR: covData.logRR[i],
        rr: Math.exp(covData.logRR[i]),
        se: Math.sqrt(covData.varLogRR[i]),
        size: 1 / covData.varLogRR[i]  // For bubble size
      });
    });
  });

  return {
    points,
    curve: modelResult.curve || [],
    referenceDose: 0,
    method: modelResult.method
  };
}

export default {
  calculateDoseCovariance,
  linearDoseResponse,
  quadraticDoseResponse,
  splineDoseResponse,
  compareDoseResponseModels,
  getDoseResponsePlotData
};
