/**
 * Heterogeneity visualization renderers
 */

import { renderPlot, scatterTrace, lineTrace, COLORS } from './plotly-utils.js';
import { getElement } from '../app/dom.js';
import { tauDL } from '../engine/tau.js';
import { poolInverseVariance } from '../engine/pooling.js';

/**
 * Render Baujat plot
 * @param {Object} result - Pooled result
 */
export function renderBaujat(result) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 2) {
    const container = getElement('baujatPlot');
    if (container) {
      container.innerHTML = '<p class="plot-message">Need ≥2 studies for Baujat plot</p>';
    }
    return;
  }

  // Color studies that contribute disproportionately to Q and influence
  const colors = active.map(s =>
    s.q > 2 && s.imp > 0.1 ? COLORS.accent : COLORS.primary
  );

  const traces = [
    {
      x: active.map(s => s.q),
      y: active.map(s => s.imp),
      mode: 'markers+text',
      type: 'scatter',
      marker: { color: colors, size: 11 },
      text: active.map(s => s.id),
      textposition: 'top center',
      textfont: { size: 9 },
      hovertemplate: '%{text}<br>Q contrib: %{x:.2f}<br>Influence: %{y:.4f}<extra></extra>'
    }
  ];

  renderPlot('baujat-plot', traces, {
    xaxis: { title: 'Contribution to Q' },
    yaxis: { title: 'Influence on Pooled Estimate' }
  });
}

/**
 * Render L'Abbé plot (for binary outcomes)
 * @param {Object} result - Pooled result
 */
export function renderLabbe(result) {
  const active = result.studies.filter(s => !s.excluded && s.raw?.p1 !== undefined);

  if (active.length < 2) {
    const container = getElement('labbePlot');
    if (container) {
      container.innerHTML = '<p class="plot-message">L\'Abbé plot requires binary outcome data</p>';
    }
    return;
  }

  // Calculate bubble sizes based on sample size
  const sizes = active.map(s => Math.sqrt((s.raw.n1 + s.raw.n2) / 10) + 5);

  const traces = [
    // Line of no effect (y = x)
    lineTrace(
      [0, 1],
      [0, 1],
      {
        line: { color: COLORS.muted, dash: 'dash', width: 1 },
        name: 'No Effect',
        hoverinfo: 'skip'
      }
    ),
    // Study bubbles
    {
      x: active.map(s => s.raw.p2),
      y: active.map(s => s.raw.p1),
      mode: 'markers+text',
      type: 'scatter',
      marker: {
        color: COLORS.primary,
        size: sizes,
        opacity: 0.7
      },
      text: active.map(s => s.id),
      textposition: 'top center',
      textfont: { size: 9 },
      name: 'Studies',
      hovertemplate: '%{text}<br>Control: %{x:.2%}<br>Treatment: %{y:.2%}<extra></extra>'
    }
  ];

  renderPlot('labbe-plot', traces, {
    xaxis: { title: 'Control Event Rate', range: [0, 1] },
    yaxis: { title: 'Treatment Event Rate', range: [0, 1] }
  });
}

/**
 * Render Radial (Galbraith) plot
 * @param {Object} result - Pooled result
 */
export function renderRadial(result) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) {
    const container = getElement('radialPlot');
    if (container) {
      container.innerHTML = '<p class="plot-message">Need ≥3 studies for radial plot</p>';
    }
    return;
  }

  const x = active.map(s => 1 / s.se); // Precision
  const y = active.map(s => s.es / s.se); // Standardized effect

  // Calculate regression line through origin
  const slope = y.reduce((acc, yi, i) => acc + yi * x[i], 0) /
                x.reduce((acc, xi) => acc + xi * xi, 0);

  const maxX = Math.max(...x) * 1.1;

  const traces = [
    // Regression line
    lineTrace(
      [0, maxX],
      [0, slope * maxX],
      {
        line: { color: COLORS.accent, width: 2 },
        name: 'Regression',
        hoverinfo: 'skip'
      }
    ),
    // ±1.96 reference lines (approximate 95% CI)
    lineTrace(
      [0, maxX],
      [1.96, 1.96],
      {
        line: { color: COLORS.muted, dash: 'dot', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),
    lineTrace(
      [0, maxX],
      [-1.96, -1.96],
      {
        line: { color: COLORS.muted, dash: 'dot', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),
    // Study points
    {
      x,
      y,
      mode: 'markers+text',
      type: 'scatter',
      marker: { color: COLORS.primary, size: 9 },
      text: active.map(s => s.id),
      textposition: 'top center',
      textfont: { size: 9 },
      name: 'Studies',
      hovertemplate: '%{text}<br>Precision: %{x:.2f}<br>Std Effect: %{y:.2f}<extra></extra>'
    }
  ];

  renderPlot('radial-plot', traces, {
    xaxis: { title: '1 / SE (Precision)' },
    yaxis: { title: 'Effect / SE (Standardized Effect)' }
  });
}

/**
 * Run and render GOSH analysis
 * @param {Object} result - Pooled result
 */
export async function runGOSHAnalysis(result) {
  if (!result) {
    alert('Please run the main analysis first');
    return;
  }

  const active = result.studies.filter(s => !s.excluded);

  if (active.length < 3 || active.length > 15) {
    alert('GOSH requires 3-15 studies (computing 2^k combinations)');
    return;
  }

  const k = active.length;
  const nComb = Math.pow(2, k) - 1;
  const results = [];

  const plotContainer = getElement('goshPlot');
  if (plotContainer) {
    plotContainer.innerHTML = '<p class="plot-message">Computing GOSH analysis...</p>';
  }

  // Generate all possible subsets
  for (let i = 1; i <= nComb; i++) {
    const subset = active.filter((_, j) => (i >> j) & 1);
    if (subset.length < 2) continue;

    const tau2 = tauDL(subset);
    const pooled = poolInverseVariance(subset.map(s => ({ ...s })), tau2, { useHKSJ: false });

    if (pooled) {
      results.push({
        es: pooled.es,
        I2: pooled.I2,
        k: subset.length
      });
    }

    // Yield to event loop periodically
    if (i % 100 === 0) {
      await new Promise(r => setTimeout(r, 0));
    }
  }

  // Render GOSH plot
  const traces = [
    {
      x: results.map(r => r.es),
      y: results.map(r => r.I2),
      mode: 'markers',
      type: 'scatter',
      marker: {
        color: results.map(r => r.k),
        colorscale: 'Viridis',
        size: 5,
        opacity: 0.6,
        colorbar: { title: 'k' }
      },
      text: results.map(r => `k=${r.k}`),
      hovertemplate: 'ES: %{x:.3f}<br>I²: %{y:.1f}%<br>%{text}<extra></extra>'
    }
  ];

  renderPlot('gosh-plot', traces, {
    xaxis: { title: 'Pooled Effect' },
    yaxis: { title: 'I² (%)', range: [0, 100] }
  });
}

/**
 * Render cumulative meta-analysis plot
 * @param {Object} result - Pooled result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderCumulative(result, isRatio) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 2) {
    const container = getElement('cumulativePlot');
    if (container) {
      container.innerHTML = '<p class="plot-message">Need ≥2 studies for cumulative analysis</p>';
    }
    return;
  }

  // Sort by year
  const sorted = [...active].sort((a, b) => {
    const yearA = parseInt(a.year) || 9999;
    const yearB = parseInt(b.year) || 9999;
    return yearA - yearB;
  });

  const cumulativeResults = [];

  for (let i = 1; i <= sorted.length; i++) {
    const subset = sorted.slice(0, i);
    const tau2 = tauDL(subset);
    const pooled = poolInverseVariance(subset.map(s => ({ ...s })), tau2, { useHKSJ: false });

    if (pooled) {
      const transform = v => isRatio ? Math.exp(v) : v;
      cumulativeResults.push({
        label: sorted[i - 1].id,
        k: i,
        es: transform(pooled.es),
        lo: transform(pooled.ciLo),
        hi: transform(pooled.ciHi)
      });
    }
  }

  const nullValue = isRatio ? 1 : 0;

  const traces = [
    // Error bars
    {
      x: cumulativeResults.map(r => r.es),
      y: cumulativeResults.map(r => r.label),
      error_x: {
        type: 'data',
        symmetric: false,
        array: cumulativeResults.map(r => r.hi - r.es),
        arrayminus: cumulativeResults.map(r => r.es - r.lo)
      },
      mode: 'markers',
      type: 'scatter',
      marker: { color: COLORS.primary, size: 9 },
      hovertemplate: '%{y}<br>ES: %{x:.3f}<extra></extra>'
    },
    // Null effect line
    lineTrace(
      [nullValue, nullValue],
      [cumulativeResults[0].label, cumulativeResults[cumulativeResults.length - 1].label],
      {
        line: { color: COLORS.accent, dash: 'dash', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    )
  ];

  renderPlot('cumulative-plot', traces, {
    xaxis: { title: 'Cumulative Effect' },
    yaxis: { autorange: 'reversed' },
    showlegend: false,
    margin: { l: 120 }
  });
}

/**
 * Render influence diagnostics
 * @param {Object} result - Pooled result
 */
export function renderInfluence(result) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) return;

  // Calculate influence diagnostics
  const stats = active.map(s => {
    const resid = (s.es - result.es) / Math.sqrt(s.vi + result.tau2);
    const dffits = resid * Math.sqrt(s.w / 100);
    const leverage = s.w / 100;
    const cookD = leverage < 1 ? resid ** 2 * leverage / (1 - leverage) : 0;

    return { id: s.id, resid, dffits, cookD, q: s.q, w: s.w };
  });

  // Update stats display
  const statsContainer = getElement('influenceStats');
  if (statsContainer) {
    statsContainer.innerHTML = `
      <div class="stat-card">
        <div class="stat-label">Max |Standardized Resid|</div>
        <div class="stat-value">${Math.max(...stats.map(s => Math.abs(s.resid))).toFixed(2)}</div>
      </div>
      <div class="stat-card">
        <div class="stat-label">Max |DFFITS|</div>
        <div class="stat-value">${Math.max(...stats.map(s => Math.abs(s.dffits))).toFixed(2)}</div>
      </div>
      <div class="stat-card">
        <div class="stat-label">Max Cook's D</div>
        <div class="stat-value">${Math.max(...stats.map(s => s.cookD)).toFixed(3)}</div>
      </div>
    `;
  }

  // Render Cook's D bar chart
  const colors = stats.map(s => s.cookD > 0.5 ? COLORS.accent : COLORS.primary);

  renderPlot('influence-plot', [{
    x: stats.map(s => s.id),
    y: stats.map(s => s.cookD),
    type: 'bar',
    marker: { color: colors },
    hovertemplate: '%{x}<br>Cook\'s D: %{y:.3f}<extra></extra>'
  }], {
    xaxis: { title: 'Study' },
    yaxis: { title: "Cook's Distance" },
    margin: { b: 80 }
  });
}
