/**
 * Plotly utilities and configuration
 */

import Plotly from 'plotly.js-dist-min';
import { PLOTLY_LAYOUT_DEFAULTS, COLORS } from '../constants/index.js';

/**
 * Track which plots have been initialized
 */
const initializedPlots = new Set();

/**
 * Render or update a Plotly plot
 * Uses Plotly.react for updates (more efficient than newPlot)
 * @param {string} elementId - Target element ID
 * @param {Array} data - Plotly data traces
 * @param {Object} layout - Plotly layout
 * @param {Object} config - Plotly config
 */
export function renderPlot(elementId, data, layout = {}, config = {}) {
  const element = document.getElementById(elementId);
  if (!element) {
    console.warn(`Plot element not found: ${elementId}`);
    return;
  }

  const fullLayout = {
    ...PLOTLY_LAYOUT_DEFAULTS,
    ...layout
  };

  const fullConfig = {
    displayModeBar: false,
    responsive: true,
    ...config
  };

  // Use Plotly.react for updates (handles both create and update efficiently)
  Plotly.react(element, data, fullLayout, fullConfig)
    .then(() => {
      initializedPlots.add(elementId);
    })
    .catch(err => {
      console.error(`Plotly render error for ${elementId}:`, err);
    });
}

/**
 * Clear a plot
 * @param {string} elementId - Target element ID
 */
export function clearPlot(elementId) {
  const element = document.getElementById(elementId);
  if (element && initializedPlots.has(elementId)) {
    Plotly.purge(element);
    initializedPlots.delete(elementId);
  }
}

/**
 * Create a scatter trace
 * @param {Array} x - X values
 * @param {Array} y - Y values
 * @param {Object} options - Trace options
 * @returns {Object} Plotly trace
 */
export function scatterTrace(x, y, options = {}) {
  return {
    x,
    y,
    mode: 'markers',
    type: 'scatter',
    marker: {
      color: COLORS.primary,
      size: 9,
      ...options.marker
    },
    ...options
  };
}

/**
 * Create a line trace
 * @param {Array} x - X values
 * @param {Array} y - Y values
 * @param {Object} options - Trace options
 * @returns {Object} Plotly trace
 */
export function lineTrace(x, y, options = {}) {
  return {
    x,
    y,
    mode: 'lines',
    type: 'scatter',
    line: {
      color: COLORS.primary,
      width: 2,
      ...options.line
    },
    ...options
  };
}

/**
 * Create a bar trace
 * @param {Array} x - X values (categories)
 * @param {Array} y - Y values
 * @param {Object} options - Trace options
 * @returns {Object} Plotly trace
 */
export function barTrace(x, y, options = {}) {
  return {
    x,
    y,
    type: 'bar',
    marker: {
      color: COLORS.primary,
      ...options.marker
    },
    ...options
  };
}

/**
 * Create error bar configuration
 * @param {Array} upper - Upper error values
 * @param {Array} lower - Lower error values (optional, defaults to symmetric)
 * @returns {Object} Error bar config
 */
export function errorBars(upper, lower = null) {
  if (lower === null) {
    return {
      type: 'data',
      array: upper,
      visible: true
    };
  }

  return {
    type: 'data',
    symmetric: false,
    array: upper,
    arrayminus: lower,
    visible: true
  };
}

/**
 * Standard x-axis configuration
 * @param {string} title - Axis title
 * @param {Object} options - Additional options
 * @returns {Object} Axis config
 */
export function xAxis(title, options = {}) {
  return {
    title,
    ...options
  };
}

/**
 * Standard y-axis configuration
 * @param {string} title - Axis title
 * @param {Object} options - Additional options
 * @returns {Object} Axis config
 */
export function yAxis(title, options = {}) {
  return {
    title,
    ...options
  };
}

/**
 * Get color based on value threshold
 * @param {number} value - Value to check
 * @param {number} threshold - Threshold value
 * @param {string} highColor - Color when above threshold
 * @param {string} lowColor - Color when below threshold
 * @returns {string} Color
 */
export function thresholdColor(value, threshold, highColor = COLORS.accent, lowColor = COLORS.primary) {
  return value > threshold ? highColor : lowColor;
}

/**
 * Generate colors for an array based on threshold
 * @param {Array} values - Array of values
 * @param {number} threshold - Threshold value
 * @returns {Array} Array of colors
 */
export function thresholdColors(values, threshold) {
  return values.map(v => thresholdColor(v, threshold));
}

/**
 * Create annotation
 * @param {number} x - X position
 * @param {number} y - Y position
 * @param {string} text - Annotation text
 * @param {Object} options - Additional options
 * @returns {Object} Annotation config
 */
export function annotation(x, y, text, options = {}) {
  return {
    x,
    y,
    text,
    showarrow: false,
    font: { size: 10 },
    ...options
  };
}
