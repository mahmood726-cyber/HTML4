/**
 * Cross-Validation and Overfitting Detection for Meta-Analysis
 *
 * Methods for model validation:
 * - Leave-one-out cross-validation (LOOCV)
 * - K-fold cross-validation
 * - Influence diagnostics
 * - Overfitting detection for meta-regression
 * - Prediction intervals
 * - External validation metrics
 *
 * Based on metaoverfit methodology and GRADE guidance
 */

import { jStat } from 'jstat';

/**
 * Leave-One-Out Cross-Validation for Meta-Analysis
 * Assesses stability and identifies influential studies
 *
 * @param {Array} studies - Studies with es and vi
 * @param {Function} poolingFn - Function to pool studies (receives array, returns {es, vi})
 * @param {Object} options - Options
 * @returns {Object} LOOCV results
 */
export function leaveOneOutCV(studies, poolingFn, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for LOOCV' };
  }

  // Full model estimate
  const fullResult = poolingFn(validStudies);
  const fullES = fullResult.es;
  const fullSE = Math.sqrt(fullResult.vi || 1 / validStudies.reduce((sum, s) =>
    sum + 1 / s.vi, 0));

  // Leave-one-out estimates
  const looResults = [];

  for (let i = 0; i < k; i++) {
    const subset = validStudies.filter((_, j) => j !== i);
    const looResult = poolingFn(subset);

    const looES = looResult.es;
    const looSE = Math.sqrt(looResult.vi || 1 / subset.reduce((sum, s) =>
      sum + 1 / s.vi, 0));

    // Influence measures
    const influence = fullES - looES;
    const relativeInfluence = fullES !== 0 ? influence / Math.abs(fullES) : influence;

    looResults.push({
      id: validStudies[i].id,
      omittedES: validStudies[i].es,
      looES,
      looSE,
      looCI: {
        lo: looES - 1.96 * looSE,
        hi: looES + 1.96 * looSE
      },
      influence,
      relativeInfluence: relativeInfluence * 100,
      influential: Math.abs(relativeInfluence) > 0.1  // >10% change
    });
  }

  // Summary statistics
  const influences = looResults.map(r => r.influence);
  const meanInfluence = influences.reduce((a, b) => a + b, 0) / k;
  const sdInfluence = Math.sqrt(
    influences.reduce((sum, inf) => sum + (inf - meanInfluence) ** 2, 0) / (k - 1)
  );

  const influentialStudies = looResults.filter(r => r.influential);

  // Stability index (coefficient of variation of LOO estimates)
  const looESValues = looResults.map(r => r.looES);
  const meanLooES = looESValues.reduce((a, b) => a + b, 0) / k;
  const sdLooES = Math.sqrt(
    looESValues.reduce((sum, es) => sum + (es - meanLooES) ** 2, 0) / (k - 1)
  );
  const stabilityCV = meanLooES !== 0 ? sdLooES / Math.abs(meanLooES) : sdLooES;

  return {
    method: 'Leave-One-Out Cross-Validation',
    k,
    fullModel: {
      es: fullES,
      se: fullSE,
      ci: {
        lo: fullES - 1.96 * fullSE,
        hi: fullES + 1.96 * fullSE
      }
    },
    results: looResults,
    summary: {
      meanInfluence,
      sdInfluence,
      nInfluential: influentialStudies.length,
      influentialStudyIds: influentialStudies.map(s => s.id),
      stabilityCV,
      isStable: stabilityCV < 0.1
    },
    interpretation: generateLOOInterpretation(influentialStudies.length, k, stabilityCV)
  };
}

function generateLOOInterpretation(nInfluential, k, cv) {
  const parts = [];

  if (nInfluential === 0) {
    parts.push('No individual study has undue influence on the pooled estimate.');
  } else if (nInfluential === 1) {
    parts.push('One study appears influential. Consider sensitivity analysis excluding it.');
  } else {
    parts.push(`${nInfluential} studies appear influential. Results may be sensitive to study inclusion.`);
  }

  if (cv < 0.05) {
    parts.push('Pooled estimate is highly stable across leave-one-out analyses.');
  } else if (cv < 0.1) {
    parts.push('Pooled estimate shows acceptable stability.');
  } else {
    parts.push('Pooled estimate shows notable variability. Interpret with caution.');
  }

  return parts.join(' ');
}

/**
 * K-Fold Cross-Validation for Meta-Regression
 * Assesses predictive performance and overfitting
 *
 * @param {Array} studies - Studies with es, vi, and covariates
 * @param {Function} regressionFn - Meta-regression function
 * @param {Array} covariateNames - Names of covariates
 * @param {Object} options - Options including k (folds)
 * @returns {Object} K-fold CV results
 */
export function kFoldCV(studies, regressionFn, covariateNames, options = {}) {
  const { nFolds = 5, nIterations = 10 } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.es != null && s.vi != null
  );
  const n = validStudies.length;

  if (n < nFolds) {
    return { error: `At least ${nFolds} studies required for ${nFolds}-fold CV` };
  }

  const allPredictions = [];
  const foldMetrics = [];

  for (let iter = 0; iter < nIterations; iter++) {
    // Shuffle studies
    const shuffled = [...validStudies].sort(() => Math.random() - 0.5);
    const foldSize = Math.floor(n / nFolds);

    for (let fold = 0; fold < nFolds; fold++) {
      // Split into train and test
      const testStart = fold * foldSize;
      const testEnd = fold === nFolds - 1 ? n : testStart + foldSize;

      const testSet = shuffled.slice(testStart, testEnd);
      const trainSet = [
        ...shuffled.slice(0, testStart),
        ...shuffled.slice(testEnd)
      ];

      // Fit model on training set
      const covariates = covariateNames.map(name =>
        trainSet.map(s => s[name])
      );

      const model = regressionFn(trainSet, covariates, covariateNames);

      if (model.error) continue;

      // Predict on test set
      for (const testStudy of testSet) {
        let predicted = model.intercept?.estimate || 0;

        if (model.slope) {
          // Univariate
          predicted += model.slope.estimate * (testStudy[covariateNames[0]] || 0);
        } else if (model.coefficients) {
          // Multivariate
          for (let c = 1; c < model.coefficients.length; c++) {
            const covName = model.coefficients[c].name;
            predicted += model.coefficients[c].estimate * (testStudy[covName] || 0);
          }
        }

        const error = testStudy.es - predicted;
        const weightedError = error * Math.sqrt(1 / testStudy.vi);

        allPredictions.push({
          iteration: iter,
          fold,
          studyId: testStudy.id,
          actual: testStudy.es,
          predicted,
          error,
          weightedError
        });
      }
    }
  }

  if (allPredictions.length === 0) {
    return { error: 'No valid predictions generated' };
  }

  // Calculate CV metrics
  const errors = allPredictions.map(p => p.error);
  const absErrors = errors.map(e => Math.abs(e));
  const squaredErrors = errors.map(e => e * e);

  const RMSE = Math.sqrt(squaredErrors.reduce((a, b) => a + b, 0) / allPredictions.length);
  const MAE = absErrors.reduce((a, b) => a + b, 0) / allPredictions.length;
  const bias = errors.reduce((a, b) => a + b, 0) / allPredictions.length;

  // Compare to null model (intercept only)
  const meanES = validStudies.reduce((sum, s) => sum + s.es, 0) / n;
  const totalSS = validStudies.reduce((sum, s) => sum + (s.es - meanES) ** 2, 0);
  const residualSS = allPredictions.reduce((sum, p) => sum + p.error ** 2, 0) / nIterations;

  const R2CV = 1 - residualSS / totalSS;

  // Shrinkage estimate
  // Compare in-sample R² to CV R²
  const inSampleSS = validStudies.reduce((sum, s) => {
    const model = regressionFn(
      validStudies,
      covariateNames.map(name => validStudies.map(st => st[name])),
      covariateNames
    );
    if (model.error) return sum;

    let predicted = model.intercept?.estimate || 0;
    if (model.slope) {
      predicted += model.slope.estimate * (s[covariateNames[0]] || 0);
    }
    return sum + (s.es - predicted) ** 2;
  }, 0);

  const R2InSample = 1 - inSampleSS / (n * totalSS / validStudies.length);
  const shrinkage = R2InSample - R2CV;

  return {
    method: 'K-Fold Cross-Validation',
    nFolds,
    nIterations,
    nStudies: n,
    metrics: {
      RMSE,
      MAE,
      bias,
      R2CV: Math.max(0, R2CV) * 100,
      R2InSample: Math.max(0, R2InSample) * 100,
      shrinkage: shrinkage * 100
    },
    overfitting: {
      detected: shrinkage > 0.1,
      severity: shrinkage > 0.2 ? 'severe' : shrinkage > 0.1 ? 'moderate' : 'minimal',
      shrinkagePct: shrinkage * 100
    },
    interpretation: generateCVInterpretation(RMSE, R2CV, shrinkage)
  };
}

function generateCVInterpretation(rmse, r2cv, shrinkage) {
  const parts = [];

  if (r2cv > 0.5) {
    parts.push('Model shows good predictive performance.');
  } else if (r2cv > 0.25) {
    parts.push('Model shows moderate predictive ability.');
  } else {
    parts.push('Model has limited predictive value.');
  }

  if (shrinkage > 0.2) {
    parts.push('Substantial overfitting detected. Consider simpler model or more data.');
  } else if (shrinkage > 0.1) {
    parts.push('Moderate overfitting. Penalized methods may help.');
  } else {
    parts.push('Minimal overfitting. Model generalizes well.');
  }

  return parts.join(' ');
}

/**
 * Influence Diagnostics for Meta-Analysis
 * Cook's distance, DFBETAS, leverage
 *
 * @param {Array} studies - Studies with es and vi
 * @param {Object} options - Options
 * @returns {Object} Influence diagnostics
 */
export function influenceDiagnostics(studies, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required' };
  }

  // Calculate weights
  const w = validStudies.map(s => 1 / s.vi);
  const sumW = w.reduce((a, b) => a + b, 0);

  // Pooled estimate
  const pooledES = w.reduce((sum, wi, i) => sum + wi * validStudies[i].es, 0) / sumW;

  // Residuals
  const residuals = validStudies.map(s => s.es - pooledES);

  // Standardized residuals
  const MSE = w.reduce((sum, wi, i) => sum + wi * residuals[i] ** 2, 0) / (k - 1);
  const standardizedResid = residuals.map((r, i) =>
    r / Math.sqrt(MSE * (1 - w[i] / sumW))
  );

  // Leverage (hat values)
  const leverage = w.map(wi => wi / sumW);

  // Cook's distance
  const cookD = validStudies.map((_, i) => {
    const h = leverage[i];
    const r = standardizedResid[i];
    return r * r * h / (1 * (1 - h) ** 2);
  });

  // DFFITS
  const dffits = validStudies.map((_, i) => {
    const h = leverage[i];
    const r = standardizedResid[i];
    return r * Math.sqrt(h / (1 - h));
  });

  // DFBETAS (influence on pooled estimate)
  const dfbetas = validStudies.map((s, i) => {
    const looSubset = validStudies.filter((_, j) => j !== i);
    const looW = looSubset.map(st => 1 / st.vi);
    const sumLooW = looW.reduce((a, b) => a + b, 0);
    const looES = looW.reduce((sum, wi, j) => sum + wi * looSubset[j].es, 0) / sumLooW;
    return (pooledES - looES) / Math.sqrt(MSE / sumW);
  });

  // Thresholds for outliers
  const cookThreshold = 4 / k;
  const dffitsThreshold = 2 * Math.sqrt(1 / k);
  const dfbetasThreshold = 2 / Math.sqrt(k);

  const diagnostics = validStudies.map((s, i) => ({
    id: s.id,
    es: s.es,
    residual: residuals[i],
    standardizedResidual: standardizedResid[i],
    leverage: leverage[i],
    cookD: cookD[i],
    dffits: dffits[i],
    dfbetas: dfbetas[i],
    outlier: Math.abs(standardizedResid[i]) > 2,
    influential: cookD[i] > cookThreshold ||
      Math.abs(dffits[i]) > dffitsThreshold ||
      Math.abs(dfbetas[i]) > dfbetasThreshold
  }));

  const outliers = diagnostics.filter(d => d.outlier);
  const influential = diagnostics.filter(d => d.influential);

  return {
    method: 'Influence Diagnostics',
    k,
    pooledES,
    thresholds: {
      cookD: cookThreshold,
      dffits: dffitsThreshold,
      dfbetas: dfbetasThreshold
    },
    diagnostics,
    summary: {
      nOutliers: outliers.length,
      outlierIds: outliers.map(d => d.id),
      nInfluential: influential.length,
      influentialIds: influential.map(d => d.id)
    },
    interpretation: generateInfluenceInterpretation(outliers.length, influential.length, k)
  };
}

function generateInfluenceInterpretation(nOutliers, nInfluential, k) {
  const parts = [];

  if (nOutliers === 0) {
    parts.push('No outlier studies detected based on standardized residuals.');
  } else {
    parts.push(`${nOutliers} potential outlier(s) identified.`);
  }

  if (nInfluential === 0) {
    parts.push('No unduly influential studies detected.');
  } else {
    const pct = (nInfluential / k * 100).toFixed(0);
    parts.push(`${nInfluential} influential study/studies (${pct}% of total).`);
  }

  return parts.join(' ');
}

/**
 * Calculate Prediction Interval
 * For estimating effect in a new study
 *
 * @param {Array} studies - Studies with es and vi
 * @param {Object} options - Options including tau2 if available
 * @returns {Object} Prediction interval
 */
export function predictionInterval(studies, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for prediction interval' };
  }

  // Calculate weights
  const w = validStudies.map(s => 1 / s.vi);
  const sumW = w.reduce((a, b) => a + b, 0);

  // Pooled estimate
  const pooledES = w.reduce((sum, wi, i) => sum + wi * validStudies[i].es, 0) / sumW;
  const pooledVar = 1 / sumW;

  // Calculate tau² if not provided
  let tau2 = options.tau2;
  if (tau2 == null) {
    const Q = w.reduce((sum, wi, i) =>
      sum + wi * (validStudies[i].es - pooledES) ** 2, 0);
    const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;
    tau2 = Math.max(0, (Q - (k - 1)) / C);
  }

  // Prediction interval variance includes:
  // 1. Variance of pooled estimate
  // 2. Between-study variance (tau²)
  // 3. Typical within-study variance (optional)
  const typicalVar = validStudies.reduce((sum, s) => sum + s.vi, 0) / k;

  const predictionVar = pooledVar + tau2;

  // Use t-distribution for small k
  const df = k - 2;
  const tCrit = df > 0 ? jStat.studentt.inv(0.975, df) : 1.96;

  const predSE = Math.sqrt(predictionVar);
  const predLo = pooledES - tCrit * predSE;
  const predHi = pooledES + tCrit * predSE;

  return {
    method: 'Prediction Interval',
    k,
    pooledES,
    tau2,
    confidenceInterval: {
      lo: pooledES - 1.96 * Math.sqrt(pooledVar),
      hi: pooledES + 1.96 * Math.sqrt(pooledVar)
    },
    predictionInterval: {
      lo: predLo,
      hi: predHi
    },
    interpretation: generatePredictionInterpretation(pooledES, predLo, predHi)
  };
}

function generatePredictionInterpretation(pooled, lo, hi) {
  const width = hi - lo;

  // Check if interval includes null
  const includesNull = (lo < 0 && hi > 0);

  let interpretation = `The effect in a new study is predicted to fall between ${lo.toFixed(2)} and ${hi.toFixed(2)}. `;

  if (includesNull) {
    interpretation += 'Note: The prediction interval includes null, suggesting heterogeneity may result in varying effects across settings.';
  } else {
    interpretation += 'The prediction interval excludes null, suggesting effects are likely to be consistently in the same direction.';
  }

  return interpretation;
}

/**
 * Overfitting Risk Assessment for Meta-Regression
 *
 * @param {Array} studies - Studies
 * @param {number} nCovariates - Number of covariates in model
 * @param {Object} options - Options
 * @returns {Object} Overfitting risk assessment
 */
export function assessOverfittingRisk(studies, nCovariates, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null);
  const k = validStudies.length;

  // Rule of thumb: need ~10 studies per covariate
  const studiesPerCovariate = k / nCovariates;
  const minRecommended = 10;

  // Events per variable (EPV) equivalent
  const epv = studiesPerCovariate;

  // Degrees of freedom
  const df = k - nCovariates - 1;

  // Risk assessment
  let riskLevel, recommendation;

  if (studiesPerCovariate < 5) {
    riskLevel = 'very high';
    recommendation = 'Severe overfitting likely. Reduce number of covariates or combine into composite.';
  } else if (studiesPerCovariate < 10) {
    riskLevel = 'high';
    recommendation = 'High overfitting risk. Consider penalized methods or variable selection.';
  } else if (studiesPerCovariate < 15) {
    riskLevel = 'moderate';
    recommendation = 'Moderate risk. Use cross-validation to assess actual overfitting.';
  } else {
    riskLevel = 'low';
    recommendation = 'Adequate sample size for number of covariates.';
  }

  // Shrinkage factor estimate (optimism-corrected)
  // Simplified van Houwelingen-le Cessie approximation
  const shrinkageFactor = Math.max(0, (df - 2) / df);

  return {
    k,
    nCovariates,
    studiesPerCovariate,
    minRecommended,
    df,
    epvEquivalent: epv,
    riskLevel,
    shrinkageFactor,
    recommendation,
    suggestions: generateOverfittingSuggestions(studiesPerCovariate, nCovariates)
  };
}

function generateOverfittingSuggestions(spc, nCov) {
  const suggestions = [];

  if (spc < 10) {
    suggestions.push('Reduce model complexity');

    if (nCov > 3) {
      suggestions.push('Limit to 2-3 most important covariates');
    }

    suggestions.push('Consider ridge/LASSO penalization');
    suggestions.push('Use cross-validation for model selection');
    suggestions.push('Report bootstrap-corrected estimates');
  }

  if (spc < 5) {
    suggestions.push('Consider pooling without covariates');
    suggestions.push('Use descriptive subgroup analysis instead');
  }

  return suggestions;
}

/**
 * Bootstrap Validation for Meta-Analysis
 * Provides optimism-corrected performance estimates
 *
 * @param {Array} studies - Studies
 * @param {Function} poolingFn - Pooling function
 * @param {Object} options - Options
 * @returns {Object} Bootstrap validation results
 */
export function bootstrapValidation(studies, poolingFn, options = {}) {
  const { nBoot = 200 } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const k = validStudies.length;

  if (k < 5) {
    return { error: 'At least 5 studies required for bootstrap validation' };
  }

  // Original estimate
  const original = poolingFn(validStudies);
  const originalES = original.es;

  // Bootstrap samples
  const bootEstimates = [];

  for (let b = 0; b < nBoot; b++) {
    // Sample with replacement
    const bootSample = [];
    for (let i = 0; i < k; i++) {
      const idx = Math.floor(Math.random() * k);
      bootSample.push(validStudies[idx]);
    }

    const bootResult = poolingFn(bootSample);
    if (bootResult && bootResult.es != null) {
      bootEstimates.push(bootResult.es);
    }
  }

  if (bootEstimates.length < nBoot * 0.9) {
    return { error: 'Too many bootstrap iterations failed' };
  }

  // Bootstrap statistics
  const sortedBoot = bootEstimates.sort((a, b) => a - b);
  const bootMean = bootEstimates.reduce((a, b) => a + b, 0) / bootEstimates.length;
  const bootVar = bootEstimates.reduce((sum, es) =>
    sum + (es - bootMean) ** 2, 0) / (bootEstimates.length - 1);

  // Bias
  const bias = bootMean - originalES;

  // Percentile CI
  const ci025 = sortedBoot[Math.floor(bootEstimates.length * 0.025)];
  const ci975 = sortedBoot[Math.floor(bootEstimates.length * 0.975)];

  // Bias-corrected accelerated (BCa) CI (simplified)
  const z0 = jStat.normal.inv(
    bootEstimates.filter(es => es < originalES).length / bootEstimates.length, 0, 1
  );

  return {
    method: 'Bootstrap Validation',
    nBootstrap: bootEstimates.length,
    original: {
      es: originalES
    },
    bootstrap: {
      mean: bootMean,
      se: Math.sqrt(bootVar),
      bias,
      percentileCI: {
        lo: ci025,
        hi: ci975
      },
      biasCorrection: z0
    },
    biasCorrectedES: originalES - bias,
    interpretation: generateBootstrapInterpretation(bias, originalES, Math.sqrt(bootVar))
  };
}

function generateBootstrapInterpretation(bias, original, se) {
  const relativeBias = original !== 0 ? Math.abs(bias / original) * 100 : 0;

  if (relativeBias < 5) {
    return 'Minimal bias detected. Original estimate appears reliable.';
  } else if (relativeBias < 10) {
    return 'Moderate bias detected. Consider reporting bias-corrected estimate.';
  } else {
    return 'Substantial bias detected. Interpret original estimate with caution.';
  }
}

export default {
  leaveOneOutCV,
  kFoldCV,
  influenceDiagnostics,
  predictionInterval,
  assessOverfittingRisk,
  bootstrapValidation
};
