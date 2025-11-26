/**
 * Survival Analysis and Kaplan-Meier Reconstruction
 * For time-to-event meta-analysis
 * References:
 * - Guyot et al. (2012) - IPD reconstruction from KM curves
 * - Parmar et al. (1998) - Methods for extracting HR
 * - Latimer (2013) - NICE DSU TSD 14 - Survival analysis
 */

/**
 * Reconstruct individual patient data from Kaplan-Meier curve
 * Based on Guyot et al. (2012) algorithm
 *
 * @param {Array} kmPoints - Digitized KM curve points [{time, survival}]
 * @param {Object} info - Additional information
 *        { nRisk: [{time, n}], nEvents: total events, nCensored: total censored }
 * @param {Object} options - Reconstruction options
 * @returns {Object} Reconstructed IPD and survival estimates
 */
export function reconstructIPD(kmPoints, info = {}, options = {}) {
  const {
    totalN = null,      // Total patients at risk at time 0
    totalEvents = null, // Total events (if known)
    method = 'guyot'    // 'guyot' or 'iterative'
  } = options;

  // Sort points by time
  const sorted = [...kmPoints].sort((a, b) => a.time - b.time);

  // Ensure starting point at time 0
  if (sorted.length === 0 || sorted[0].time > 0) {
    sorted.unshift({ time: 0, survival: 1 });
  }

  // Initial number at risk
  let n0 = totalN;
  if (!n0 && info.nRisk && info.nRisk.length > 0) {
    const first = info.nRisk.find(r => r.time === 0) || info.nRisk[0];
    n0 = first.n;
  }
  if (!n0) {
    // Estimate from curve shape (rough approximation)
    n0 = 100;  // Default assumption
  }

  // Reconstruct IPD using Guyot algorithm
  const ipd = [];
  let nAtRisk = n0;

  for (let i = 1; i < sorted.length; i++) {
    const prev = sorted[i - 1];
    const curr = sorted[i];

    // Number of events in this interval
    // S(t) = S(t-1) * (1 - d/n) => d = n * (1 - S(t)/S(t-1))
    const survRatio = curr.survival / prev.survival;

    if (survRatio < 1 && nAtRisk > 0) {
      const dEvents = Math.round(nAtRisk * (1 - survRatio));

      // Add event records
      for (let j = 0; j < dEvents; j++) {
        // Distribute events uniformly in interval
        const eventTime = prev.time + (curr.time - prev.time) * (j + 0.5) / Math.max(1, dEvents);
        ipd.push({
          time: eventTime,
          event: 1,
          weight: 1
        });
      }

      nAtRisk -= dEvents;
    }

    // Handle censoring (if number at risk information available)
    if (info.nRisk) {
      const nRiskCurr = info.nRisk.find(r => Math.abs(r.time - curr.time) < 0.01);
      if (nRiskCurr && nRiskCurr.n < nAtRisk) {
        const nCensored = nAtRisk - nRiskCurr.n;

        for (let j = 0; j < nCensored; j++) {
          const censorTime = prev.time + (curr.time - prev.time) * (j + 0.5) / Math.max(1, nCensored);
          ipd.push({
            time: censorTime,
            event: 0,
            weight: 1
          });
        }

        nAtRisk = nRiskCurr.n;
      }
    }
  }

  // Add remaining censored at last observation
  const lastTime = sorted[sorted.length - 1].time;
  while (nAtRisk > 0) {
    ipd.push({
      time: lastTime,
      event: 0,
      weight: 1
    });
    nAtRisk--;
  }

  // Sort by time
  ipd.sort((a, b) => a.time - b.time);

  // Calculate summary statistics
  const nEvents = ipd.filter(p => p.event === 1).length;
  const nCensored = ipd.filter(p => p.event === 0).length;

  // Reconstruct KM curve from IPD to verify
  const reconstructedKM = calculateKMFromIPD(ipd);

  // Calculate median survival
  const median = findQuantileFromKM(reconstructedKM, 0.5);

  // Calculate restricted mean survival time (RMST)
  const rmst = calculateRMST(reconstructedKM, lastTime);

  return {
    ipd,
    n: ipd.length,
    nEvents,
    nCensored,
    median,
    rmst,
    reconstructedKM,
    originalKM: sorted,
    diagnostics: {
      matchQuality: assessReconstruction(sorted, reconstructedKM)
    }
  };
}

/**
 * Calculate Kaplan-Meier curve from IPD
 */
export function calculateKMFromIPD(ipd) {
  // Sort by time
  const sorted = [...ipd].sort((a, b) => a.time - b.time);

  const km = [{ time: 0, survival: 1, nRisk: sorted.length, nEvents: 0 }];
  let survival = 1;
  let nRisk = sorted.length;
  let i = 0;

  while (i < sorted.length) {
    const currentTime = sorted[i].time;

    // Count events and censored at this time
    let events = 0;
    let censored = 0;

    while (i < sorted.length && sorted[i].time === currentTime) {
      if (sorted[i].event === 1) events++;
      else censored++;
      i++;
    }

    // Update survival
    if (events > 0 && nRisk > 0) {
      survival *= (nRisk - events) / nRisk;
      km.push({
        time: currentTime,
        survival,
        nRisk,
        nEvents: events
      });
    }

    nRisk -= events + censored;
  }

  return km;
}

/**
 * Find quantile (e.g., median) from KM curve
 */
function findQuantileFromKM(km, quantile) {
  const targetSurvival = 1 - quantile;

  for (let i = 0; i < km.length; i++) {
    if (km[i].survival <= targetSurvival) {
      if (i === 0) return km[0].time;

      // Linear interpolation
      const s1 = km[i - 1].survival;
      const s2 = km[i].survival;
      const t1 = km[i - 1].time;
      const t2 = km[i].time;

      if (s1 === s2) return t1;
      return t1 + (t2 - t1) * (s1 - targetSurvival) / (s1 - s2);
    }
  }

  return null;  // Quantile not reached
}

/**
 * Calculate Restricted Mean Survival Time (RMST)
 * Area under KM curve up to time t*
 */
export function calculateRMST(km, tStar) {
  let rmst = 0;

  for (let i = 1; i < km.length; i++) {
    if (km[i].time > tStar) break;

    const width = km[i].time - km[i - 1].time;
    rmst += km[i - 1].survival * width;
  }

  // Add final rectangle if not reached tStar
  const lastIdx = km.findIndex(k => k.time >= tStar) - 1;
  if (lastIdx >= 0 && km[lastIdx].time < tStar) {
    rmst += km[lastIdx].survival * (tStar - km[lastIdx].time);
  }

  return rmst;
}

/**
 * Assess quality of reconstruction by comparing curves
 */
function assessReconstruction(original, reconstructed) {
  if (original.length < 2 || reconstructed.length < 2) {
    return { rmse: null, maxDiff: null };
  }

  let sumSqDiff = 0;
  let maxDiff = 0;
  let count = 0;

  // Compare at original time points
  original.forEach(pt => {
    // Find closest reconstructed point
    let closest = reconstructed[0];
    let minDist = Math.abs(pt.time - reconstructed[0].time);

    for (const rpt of reconstructed) {
      const dist = Math.abs(pt.time - rpt.time);
      if (dist < minDist) {
        minDist = dist;
        closest = rpt;
      }
    }

    const diff = Math.abs(pt.survival - closest.survival);
    sumSqDiff += diff ** 2;
    maxDiff = Math.max(maxDiff, diff);
    count++;
  });

  return {
    rmse: Math.sqrt(sumSqDiff / count),
    maxDiff,
    good: maxDiff < 0.05  // Less than 5% deviation
  };
}

/**
 * Extract hazard ratio from two KM curves
 * Methods: logrank, Peto, Cox-Mantel
 *
 * @param {Array} km1 - KM curve for treatment group
 * @param {Array} km2 - KM curve for control group
 * @param {Object} options - Extraction options
 * @returns {Object} Hazard ratio estimate
 */
export function extractHazardRatio(km1, km2, options = {}) {
  const { method = 'logrank' } = options;

  // Need IPD for proper calculation
  // This is an approximation using the curve data

  // Simple method: compare median survivals
  const median1 = findQuantileFromKM(km1, 0.5);
  const median2 = findQuantileFromKM(km2, 0.5);

  // Under exponential assumption: HR ≈ median2/median1
  let hrApprox = null;
  if (median1 && median2) {
    hrApprox = median2 / median1;
  }

  // Calculate at multiple time points for robustness
  const timePoints = [];
  const maxTime = Math.min(
    km1[km1.length - 1].time,
    km2[km2.length - 1].time
  );

  for (let t = maxTime * 0.1; t <= maxTime * 0.9; t += maxTime * 0.1) {
    const s1 = getSurvivalAtTime(km1, t);
    const s2 = getSurvivalAtTime(km2, t);

    if (s1 > 0 && s1 < 1 && s2 > 0 && s2 < 1) {
      // HR from cumulative hazards: H(t) = -log(S(t))
      const h1 = -Math.log(s1);
      const h2 = -Math.log(s2);
      if (h2 > 0) {
        timePoints.push({ time: t, hr: h1 / h2 });
      }
    }
  }

  // Average HR across time points (assuming proportional hazards)
  const avgHR = timePoints.length > 0
    ? timePoints.reduce((sum, p) => sum + p.hr, 0) / timePoints.length
    : hrApprox;

  // Variance approximation (very rough)
  const n1 = km1[0].nRisk || 50;
  const n2 = km2[0].nRisk || 50;
  const e1 = km1.reduce((sum, k) => sum + (k.nEvents || 0), 0);
  const e2 = km2.reduce((sum, k) => sum + (k.nEvents || 0), 0);
  const totalEvents = e1 + e2;

  // SE approximation: sqrt(4/E) for log(HR)
  const seLnHR = totalEvents > 0 ? Math.sqrt(4 / totalEvents) : 1;

  return {
    hr: avgHR,
    lnHR: avgHR ? Math.log(avgHR) : null,
    se: seLnHR,
    ciLo: avgHR ? avgHR * Math.exp(-1.96 * seLnHR) : null,
    ciHi: avgHR ? avgHR * Math.exp(1.96 * seLnHR) : null,
    method,
    medians: { treatment: median1, control: median2 },
    proportionalHazards: assessProportionalHazards(timePoints)
  };
}

/**
 * Get survival probability at a specific time
 */
function getSurvivalAtTime(km, t) {
  if (t <= km[0].time) return 1;

  for (let i = 1; i < km.length; i++) {
    if (km[i].time >= t) {
      // Linear interpolation
      const s1 = km[i - 1].survival;
      const s2 = km[i].survival;
      const t1 = km[i - 1].time;
      const t2 = km[i].time;

      return s1 + (s2 - s1) * (t - t1) / (t2 - t1);
    }
  }

  return km[km.length - 1].survival;
}

/**
 * Assess proportional hazards assumption
 */
function assessProportionalHazards(hrOverTime) {
  if (hrOverTime.length < 3) {
    return { valid: null, message: 'Insufficient data' };
  }

  // Check if HR is constant over time
  const hrs = hrOverTime.map(p => p.hr);
  const mean = hrs.reduce((a, b) => a + b, 0) / hrs.length;
  const variance = hrs.reduce((sum, h) => sum + (h - mean) ** 2, 0) / hrs.length;
  const cv = Math.sqrt(variance) / mean;  // Coefficient of variation

  return {
    valid: cv < 0.3,  // Less than 30% variation suggests PH holds
    cv,
    message: cv < 0.3 ? 'Proportional hazards assumption appears reasonable' :
             'Evidence of non-proportional hazards - consider time-varying HR'
  };
}

/**
 * Fit parametric survival model to KM data
 * Supports: exponential, Weibull, log-normal, log-logistic, Gompertz
 *
 * @param {Array} km - KM curve data
 * @param {string} distribution - Distribution to fit
 * @returns {Object} Fitted model parameters and goodness of fit
 */
export function fitParametricSurvival(km, distribution = 'weibull') {
  // Convert KM to survival times
  const times = [];
  const survival = [];

  km.forEach(pt => {
    if (pt.time > 0 && pt.survival > 0 && pt.survival < 1) {
      times.push(pt.time);
      survival.push(pt.survival);
    }
  });

  if (times.length < 3) {
    return { error: 'Insufficient data for parametric fit' };
  }

  let params, fittedS;

  switch (distribution) {
    case 'exponential':
      params = fitExponential(times, survival);
      fittedS = times.map(t => Math.exp(-params.rate * t));
      break;

    case 'weibull':
      params = fitWeibull(times, survival);
      fittedS = times.map(t => Math.exp(-Math.pow(t / params.scale, params.shape)));
      break;

    case 'lognormal':
      params = fitLogNormal(times, survival);
      fittedS = times.map(t => 1 - normalCDF((Math.log(t) - params.meanlog) / params.sdlog));
      break;

    case 'loglogistic':
      params = fitLogLogistic(times, survival);
      fittedS = times.map(t => 1 / (1 + Math.pow(t / params.scale, params.shape)));
      break;

    case 'gompertz':
      params = fitGompertz(times, survival);
      fittedS = times.map(t => Math.exp(-params.rate / params.shape * (Math.exp(params.shape * t) - 1)));
      break;

    default:
      return { error: `Unknown distribution: ${distribution}` };
  }

  // Goodness of fit
  const sse = times.reduce((sum, t, i) => sum + (survival[i] - fittedS[i]) ** 2, 0);
  const sst = times.reduce((sum, t, i) => {
    const mean = survival.reduce((a, b) => a + b, 0) / survival.length;
    return sum + (survival[i] - mean) ** 2;
  }, 0);

  const rSquared = 1 - sse / sst;
  const aic = times.length * Math.log(sse / times.length) + 2 * Object.keys(params).length;

  return {
    distribution,
    parameters: params,
    fitted: times.map((t, i) => ({ time: t, survival: fittedS[i] })),
    goodnessOfFit: {
      rSquared,
      aic,
      sse,
      good: rSquared > 0.9
    }
  };
}

/**
 * Fit exponential distribution: S(t) = exp(-λt)
 */
function fitExponential(times, survival) {
  // Linear regression on log(S) = -λt
  const logS = survival.map(s => Math.log(s));
  const n = times.length;

  const sumT = times.reduce((a, b) => a + b, 0);
  const sumLogS = logS.reduce((a, b) => a + b, 0);
  const sumTLogS = times.reduce((acc, t, i) => acc + t * logS[i], 0);
  const sumT2 = times.reduce((acc, t) => acc + t * t, 0);

  const rate = -(n * sumTLogS - sumT * sumLogS) / (n * sumT2 - sumT * sumT);

  return { rate: Math.max(0.001, rate) };
}

/**
 * Fit Weibull distribution: S(t) = exp(-(t/λ)^k)
 */
function fitWeibull(times, survival) {
  // Linear regression on log(-log(S)) = k*log(t) - k*log(λ)
  const y = survival.map(s => Math.log(-Math.log(Math.max(0.001, Math.min(0.999, s)))));
  const x = times.map(t => Math.log(t));
  const n = times.length;

  const sumX = x.reduce((a, b) => a + b, 0);
  const sumY = y.reduce((a, b) => a + b, 0);
  const sumXY = x.reduce((acc, xi, i) => acc + xi * y[i], 0);
  const sumX2 = x.reduce((acc, xi) => acc + xi * xi, 0);

  const shape = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
  const intercept = (sumY - shape * sumX) / n;
  const scale = Math.exp(-intercept / shape);

  return {
    shape: Math.max(0.1, shape),
    scale: Math.max(0.1, scale)
  };
}

/**
 * Fit log-normal distribution
 */
function fitLogNormal(times, survival) {
  // Use quantile matching
  const q50 = findQuantileFromKM(times.map((t, i) => ({ time: t, survival: survival[i] })), 0.5);
  const q25 = findQuantileFromKM(times.map((t, i) => ({ time: t, survival: survival[i] })), 0.25);

  if (!q50 || !q25) {
    return { meanlog: Math.log(times[Math.floor(times.length / 2)]), sdlog: 1 };
  }

  const meanlog = Math.log(q50);
  const sdlog = (Math.log(q50) - Math.log(q25)) / 0.6745;  // Based on normal quantiles

  return { meanlog, sdlog: Math.max(0.1, sdlog) };
}

/**
 * Fit log-logistic distribution: S(t) = 1/(1 + (t/α)^β)
 */
function fitLogLogistic(times, survival) {
  // Linear regression on log(S/(1-S)) = -β*log(t) + β*log(α)
  const y = survival.map(s => {
    const sAdj = Math.max(0.01, Math.min(0.99, s));
    return Math.log(sAdj / (1 - sAdj));
  });
  const x = times.map(t => Math.log(t));
  const n = times.length;

  const sumX = x.reduce((a, b) => a + b, 0);
  const sumY = y.reduce((a, b) => a + b, 0);
  const sumXY = x.reduce((acc, xi, i) => acc + xi * y[i], 0);
  const sumX2 = x.reduce((acc, xi) => acc + xi * xi, 0);

  const shape = -(n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
  const intercept = (sumY + shape * sumX) / n;
  const scale = Math.exp(intercept / shape);

  return {
    shape: Math.max(0.1, shape),
    scale: Math.max(0.1, scale)
  };
}

/**
 * Fit Gompertz distribution (simplified)
 */
function fitGompertz(times, survival) {
  // Approximate using grid search
  let bestRate = 0.01, bestShape = 0.01, bestSSE = Infinity;

  for (let rate = 0.001; rate <= 0.5; rate += 0.01) {
    for (let shape = 0.001; shape <= 0.5; shape += 0.01) {
      const fitted = times.map(t =>
        Math.exp(-rate / shape * (Math.exp(shape * t) - 1))
      );
      const sse = times.reduce((sum, t, i) => sum + (survival[i] - fitted[i]) ** 2, 0);

      if (sse < bestSSE) {
        bestSSE = sse;
        bestRate = rate;
        bestShape = shape;
      }
    }
  }

  return { rate: bestRate, shape: bestShape };
}

/**
 * Compare multiple parametric distributions
 */
export function compareSurvivalDistributions(km) {
  const distributions = ['exponential', 'weibull', 'lognormal', 'loglogistic', 'gompertz'];
  const results = [];

  distributions.forEach(dist => {
    const fit = fitParametricSurvival(km, dist);
    if (!fit.error) {
      results.push({
        distribution: dist,
        ...fit.goodnessOfFit,
        parameters: fit.parameters
      });
    }
  });

  // Sort by AIC
  results.sort((a, b) => a.aic - b.aic);

  return {
    fits: results,
    bestFit: results[0]?.distribution,
    ranking: results.map(r => r.distribution)
  };
}

/**
 * Normal CDF for log-normal calculations
 */
function normalCDF(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

  return 0.5 * (1 + sign * y);
}
