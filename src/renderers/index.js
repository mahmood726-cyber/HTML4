/**
 * Renderers Module - Main exports and orchestration
 */

// Re-export individual renderers
export { renderPlot, clearPlot } from './plotly-utils.js';
export { renderForest, renderLOO } from './forest.js';
export { renderBias, renderFunnel, renderContourFunnel, renderSunset } from './bias.js';
export {
  renderBaujat,
  renderLabbe,
  renderRadial,
  runGOSHAnalysis,
  renderCumulative,
  renderInfluence
} from './heterogeneity.js';
export { renderGRADE, renderRoB } from './grade.js';
export { renderEconomic } from './economic.js';
export { renderTSA, getTSAInterpretation } from './tsa.js';

// Import for orchestration
import { renderForest, renderLOO } from './forest.js';
import { renderBias } from './bias.js';
import { renderBaujat, renderLabbe, renderRadial, renderCumulative, renderInfluence } from './heterogeneity.js';
import { renderGRADE, renderRoB } from './grade.js';
import { renderEconomic } from './economic.js';
import { renderTSA } from './tsa.js';

/**
 * Render all visualizations
 * @param {Object} result - Pooled meta-analysis result
 * @param {Object} state - Application state
 */
export function renderAll(result, state) {
  if (!result) return;

  const { isRatio, metric, type } = result;

  // Core visualizations
  renderForest(result, isRatio);
  renderLOO(result, isRatio);

  // Heterogeneity plots
  renderBaujat(result);
  if (type === 'binary') {
    renderLabbe(result);
  }
  renderRadial(result);
  renderCumulative(result, isRatio);
  renderInfluence(result);

  // Publication bias
  renderBias(result);

  // TSA
  renderTSA(result);

  // GRADE and RoB
  renderGRADE(result);
  renderRoB(result);

  // Economic (only for binary outcomes)
  if (['OR', 'RR', 'RD'].includes(metric)) {
    renderEconomic(result, metric);
  }
}

/**
 * Clear all plots
 */
export function clearAll() {
  const plotIds = [
    'forest-plot',
    'funnel-plot',
    'contour-plot',
    'sunset-plot',
    'baujat-plot',
    'labbe-plot',
    'radial-plot',
    'gosh-plot',
    'tsa-plot',
    'loo-plot',
    'cumulative-plot',
    'influence-plot',
    'rob-plot'
  ];

  const { clearPlot } = require('./plotly-utils.js');
  plotIds.forEach(id => clearPlot(id));
}
