/**
 * Trial Sequential Analysis renderer
 */

import { renderPlot, lineTrace, COLORS } from './plotly-utils.js';
import { getElement, getText, setValue, getValue } from '../app/dom.js';
import { calculateTSA, getTSAPlotData } from '../engine/tsa.js';

/**
 * Render TSA analysis
 * @param {Object} result - Pooled result
 */
export function renderTSA(result) {
  // Get TSA parameters
  const alpha = parseFloat(getValue('tsaAlpha')) || 0.05;
  const beta = parseFloat(getValue('tsaBeta')) || 0.20;
  const rrr = (parseFloat(getValue('tsaRrr')) || 20) / 100;
  const cer = parseFloat(getValue('cer')) || 0.10;

  // Calculate TSA
  const tsaResult = calculateTSA(result.studies, { alpha, beta, rrr, cer });

  if (!tsaResult) {
    updateTSADisplay(null);
    return;
  }

  updateTSADisplay(tsaResult);
  renderTSAPlot(tsaResult);
}

/**
 * Update TSA display values
 * @param {Object|null} tsa - TSA result
 */
function updateTSADisplay(tsa) {
  const setText = (id, value) => {
    const el = document.getElementById(id);
    if (el) el.textContent = value;
  };

  if (!tsa) {
    setText('tsa-ris', '—');
    setText('tsa-accrued', '—');
    setText('tsa-pct', '—');
    setText('tsa-d2', '—');
    setText('tsa-conclusion', '—');
    setText('tsa-boundary', '—');
    return;
  }

  setText('tsa-ris', tsa.RIS.toLocaleString());
  setText('tsa-accrued', tsa.accrued.toLocaleString());
  setText('tsa-pct', `${tsa.pctRIS}% of RIS`);
  setText('tsa-d2', `${tsa.D2}%`);
  setText('tsa-conclusion', tsa.conclusion);
  setText('tsa-boundary', tsa.boundary);

  // Color the conclusion based on result
  const conclusionEl = document.getElementById('tsa-conclusion');
  if (conclusionEl) {
    if (tsa.crossedBenefit) {
      conclusionEl.style.color = 'var(--success)';
    } else if (tsa.crossedHarm) {
      conclusionEl.style.color = 'var(--warning)';
    } else if (tsa.crossedFutility) {
      conclusionEl.style.color = 'var(--ink-muted)';
    } else {
      conclusionEl.style.color = 'var(--ink)';
    }
  }
}

/**
 * Render TSA plot
 * @param {Object} tsa - TSA result
 */
function renderTSAPlot(tsa) {
  const plotData = getTSAPlotData(tsa);
  if (!plotData) return;

  const traces = [
    // Traditional significance boundary (dashed)
    lineTrace(
      [0, plotData.maxN],
      [plotData.traditionalZ, plotData.traditionalZ],
      {
        line: { color: COLORS.muted, dash: 'dash', width: 1 },
        name: 'Traditional (Z=1.96)',
        hoverinfo: 'skip'
      }
    ),
    lineTrace(
      [0, plotData.maxN],
      [-plotData.traditionalZ, -plotData.traditionalZ],
      {
        line: { color: COLORS.muted, dash: 'dash', width: 1 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),

    // Upper monitoring boundary
    lineTrace(
      plotData.cumN,
      plotData.monitoringUpper,
      {
        line: { color: COLORS.warning, width: 2 },
        name: 'Monitoring Boundary',
        hoverinfo: 'skip'
      }
    ),
    // Lower monitoring boundary
    lineTrace(
      plotData.cumN,
      plotData.monitoringLower,
      {
        line: { color: COLORS.warning, width: 2 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),

    // Upper futility boundary
    lineTrace(
      plotData.cumN,
      plotData.futilityUpper,
      {
        line: { color: COLORS.success, width: 2 },
        name: 'Futility Boundary',
        hoverinfo: 'skip'
      }
    ),
    // Lower futility boundary
    lineTrace(
      plotData.cumN,
      plotData.futilityLower,
      {
        line: { color: COLORS.success, width: 2 },
        showlegend: false,
        hoverinfo: 'skip'
      }
    ),

    // Cumulative Z-curve
    {
      x: plotData.cumN,
      y: plotData.cumZ,
      mode: 'lines+markers',
      type: 'scatter',
      line: { color: COLORS.info, width: 3 },
      marker: { size: 6 },
      name: 'Z-curve',
      hovertemplate: 'N: %{x}<br>Z: %{y:.2f}<extra></extra>'
    },

    // Required Information Size line
    lineTrace(
      [plotData.RIS, plotData.RIS],
      [-6, 6],
      {
        line: { color: COLORS.accent, dash: 'dot', width: 2 },
        name: 'RIS',
        hoverinfo: 'skip'
      }
    )
  ];

  renderPlot('tsa-plot', traces, {
    xaxis: {
      title: 'Cumulative Sample Size',
      range: [0, plotData.maxN]
    },
    yaxis: {
      title: 'Cumulative Z-score',
      range: [-6, 6]
    },
    showlegend: true,
    legend: {
      x: 0,
      y: 1,
      bgcolor: 'rgba(255,255,255,0.8)',
      font: { size: 10 }
    }
  });
}

/**
 * Get TSA interpretation text
 * @param {Object} tsa - TSA result
 * @returns {string} Interpretation text
 */
export function getTSAInterpretation(tsa) {
  if (!tsa) return 'Unable to perform TSA analysis.';

  const pctAcquired = parseFloat(tsa.pctRIS);

  let interpretation = `Required Information Size: ${tsa.RIS.toLocaleString()} participants. `;
  interpretation += `Currently accrued: ${tsa.accrued.toLocaleString()} (${tsa.pctRIS}%). `;

  if (tsa.crossedBenefit) {
    interpretation += 'The cumulative Z-curve has crossed the monitoring boundary for benefit. ';
    interpretation += 'There is sufficient evidence to conclude that the intervention is effective, ';
    interpretation += 'accounting for repeated testing and random error.';
  } else if (tsa.crossedHarm) {
    interpretation += 'The cumulative Z-curve has crossed the monitoring boundary for harm. ';
    interpretation += 'There is sufficient evidence to conclude that the intervention is harmful.';
  } else if (tsa.crossedFutility) {
    interpretation += 'The cumulative Z-curve has entered the futility region. ';
    interpretation += 'It is unlikely that a statistically significant effect will be found ';
    interpretation += 'even with additional studies.';
  } else if (pctAcquired >= 100) {
    interpretation += 'The required information size has been reached without crossing monitoring boundaries. ';
    interpretation += 'The evidence is inconclusive for the anticipated effect size.';
  } else {
    interpretation += 'The analysis is inconclusive. ';
    interpretation += `Additional studies totaling approximately ${(tsa.RIS - tsa.accrued).toLocaleString()} `;
    interpretation += 'more participants may be needed.';
  }

  return interpretation;
}
