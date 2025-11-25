/**
 * Publication bias visualization renderers
 */

import { renderPlot, scatterTrace, lineTrace, COLORS } from './plotly-utils.js';
import { getElement, setText } from '../app/dom.js';
import { formatNumber, formatPValue, formatEffectSize } from '../utils/index.js';
import {
  eggerTest,
  beggTest,
  failSafeN,
  trimAndFill,
  getFunnelPlotData
} from '../engine/bias.js';

/**
 * Render all publication bias outputs
 * @param {Object} result - Pooled result
 */
export function renderBias(result) {
  renderBiasStats(result);
  renderFunnel(result, result.isRatio);
  renderContourFunnel(result, result.isRatio);
}

/**
 * Render bias test statistics
 * @param {Object} result - Pooled result
 */
function renderBiasStats(result) {
  const egger = eggerTest(result.studies);
  const begg = beggTest(result.studies);
  const fsn = failSafeN(result.studies);
  const tf = trimAndFill(result.studies);

  setText('eggerP', formatPValue(egger.p));
  setText('eggerInt', egger.intercept !== null ? `Int: ${formatNumber(egger.intercept)}` : '');

  setText('beggP', formatPValue(begg.p));
  setText('beggTau', begg.tau !== null ? `τ = ${formatNumber(begg.tau)}` : '');

  setText('fsnVal', fsn.toString());
  setText('tfMissing', `${tf.k0} imputed`);

  if (tf.adjusted) {
    const adjES = formatEffectSize(tf.adjusted.es, result.isRatio);
    setText('tfAdjusted', `Adj: ${adjES}`);
  } else {
    setText('tfAdjusted', '—');
  }
}

/**
 * Render funnel plot
 * @param {Object} result - Pooled result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderFunnel(result, isRatio) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) {
    const container = getElement('funnelPlot');
    if (container) {
      container.innerHTML = '<p class="plot-message">Need ≥3 studies for funnel plot</p>';
    }
    return;
  }

  const plotData = getFunnelPlotData(result.studies, result, isRatio);
  if (!plotData) return;

  const transform = v => isRatio ? Math.exp(v) : v;
  const pooledES = transform(result.es);
  const maxSE = plotData.maxSE;

  const traces = [
    // 95% CI boundaries
    lineTrace(
      plotData.ciLines.lower,
      plotData.ciLines.se,
      {
        line: { color: 'rgba(13,115,119,0.3)', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),
    lineTrace(
      plotData.ciLines.upper,
      plotData.ciLines.se,
      {
        line: { color: 'rgba(13,115,119,0.3)', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),
    // Pooled effect line
    lineTrace(
      [pooledES, pooledES],
      [0, maxSE],
      {
        line: { color: COLORS.accent, width: 2, dash: 'dash' },
        name: 'Pooled'
      }
    ),
    // Study points
    scatterTrace(
      plotData.studies.map(s => s.x),
      plotData.studies.map(s => s.y),
      {
        marker: { color: COLORS.primary, size: 9 },
        text: plotData.studies.map(s => s.id),
        name: 'Studies',
        hovertemplate: '%{text}<br>Effect: %{x:.3f}<br>SE: %{y:.3f}<extra></extra>'
      }
    )
  ];

  renderPlot('funnel-plot', traces, {
    xaxis: { title: 'Effect Size' },
    yaxis: { title: 'Standard Error', autorange: 'reversed' },
    showlegend: false
  });
}

/**
 * Render contour-enhanced funnel plot
 * @param {Object} result - Pooled result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderContourFunnel(result, isRatio) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) return;

  const maxSE = Math.max(...active.map(e => e.se)) * 1.1;
  const seRange = Array.from({ length: 50 }, (_, i) => (i * maxSE / 49) + 0.01);

  const transform = v => isRatio ? Math.exp(v) : v;

  // Significance contours
  const contours = [
    { z: 2.576, label: 'p<0.01', color: 'rgba(220, 38, 38, 0.15)' },
    { z: 1.96, label: 'p<0.05', color: 'rgba(234, 179, 8, 0.15)' },
    { z: 1.645, label: 'p<0.10', color: 'rgba(34, 197, 94, 0.15)' }
  ];

  const traces = [];

  // Add contour regions
  contours.forEach(({ z, label, color }) => {
    traces.push({
      x: seRange.map(se => transform(-z * se)).concat(seRange.map(se => transform(z * se)).reverse()),
      y: seRange.concat([...seRange].reverse()),
      fill: 'toself',
      fillcolor: color,
      line: { width: 0 },
      name: label,
      hoverinfo: 'skip'
    });
  });

  // Add zero line
  traces.push(lineTrace(
    [transform(0), transform(0)],
    [0, maxSE],
    {
      line: { color: '#666', width: 1, dash: 'dot' },
      showlegend: false,
      hoverinfo: 'skip'
    }
  ));

  // Add study points
  traces.push(scatterTrace(
    active.map(e => transform(e.es)),
    active.map(e => e.se),
    {
      marker: { color: COLORS.primary, size: 9 },
      text: active.map(e => e.id),
      name: 'Studies',
      hovertemplate: '%{text}<br>Effect: %{x:.3f}<br>SE: %{y:.3f}<extra></extra>'
    }
  ));

  renderPlot('contour-plot', traces, {
    xaxis: { title: 'Effect Size' },
    yaxis: { title: 'Standard Error', autorange: 'reversed' },
    showlegend: true,
    legend: { x: 0.02, y: 0.98 }
  });
}

/**
 * Render sunset (power-enhanced) funnel plot
 * @param {Object} result - Pooled result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderSunset(result, isRatio) {
  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) return;

  // Sunset plot shows studies colored by statistical power
  const transform = v => isRatio ? Math.exp(v) : v;

  // Calculate power for each study (simplified)
  const studiesWithPower = active.map(s => {
    const z = Math.abs(s.es / s.se);
    const power = Math.min(0.99, 1 - 0.5 * Math.exp(-0.5 * z * z)); // Approximate
    return { ...s, power };
  });

  const colors = studiesWithPower.map(s => {
    if (s.power >= 0.8) return COLORS.success;
    if (s.power >= 0.5) return '#eab308';
    return COLORS.accent;
  });

  const traces = [
    scatterTrace(
      studiesWithPower.map(s => transform(s.es)),
      studiesWithPower.map(s => s.se),
      {
        marker: {
          color: colors,
          size: 10,
          line: { color: '#fff', width: 1 }
        },
        text: studiesWithPower.map(s => `${s.id} (Power: ${(s.power * 100).toFixed(0)}%)`),
        hovertemplate: '%{text}<extra></extra>'
      }
    )
  ];

  renderPlot('sunset-plot', traces, {
    xaxis: { title: 'Effect Size' },
    yaxis: { title: 'Standard Error', autorange: 'reversed' },
    showlegend: false,
    annotations: [
      { x: 0.02, y: 0.98, xref: 'paper', yref: 'paper', text: '● High power (≥80%)', showarrow: false, font: { size: 10, color: COLORS.success } },
      { x: 0.02, y: 0.92, xref: 'paper', yref: 'paper', text: '● Medium power (50-80%)', showarrow: false, font: { size: 10, color: '#eab308' } },
      { x: 0.02, y: 0.86, xref: 'paper', yref: 'paper', text: '● Low power (<50%)', showarrow: false, font: { size: 10, color: COLORS.accent } }
    ]
  });
}
