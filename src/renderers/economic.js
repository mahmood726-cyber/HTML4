/**
 * Health economic outputs renderer
 */

import { getElement, setText, setHTML, getValue } from '../app/dom.js';
import { formatNumber, clamp } from '../utils/index.js';
import { WTP_THRESHOLD } from '../constants/index.js';

/**
 * Render economic outputs
 * @param {Object} result - Pooled result
 * @param {string} metric - Effect measure
 */
export function renderEconomic(result, metric) {
  if (!['OR', 'RR', 'RD'].includes(metric)) {
    setText('nntValue', 'N/A');
    setText('nntLabel', 'Requires binary outcome');
    return;
  }

  // Get economic inputs
  const cer = parseFloat(getValue('cer')) || 0.1;
  const costEvent = parseFloat(getValue('costEvent')) || 10000;
  const costTreatment = parseFloat(getValue('costTreatment')) || 500;
  const qalyLoss = parseFloat(getValue('qalyLoss')) || 0.5;

  // Calculate EER from effect size
  const es = result.es;
  let eer;

  if (metric === 'OR') {
    const or = Math.exp(es);
    eer = (cer * or) / (1 - cer + cer * or);
  } else if (metric === 'RR') {
    eer = cer * Math.exp(es);
  } else {
    // RD
    eer = cer + es;
  }

  // Clamp to valid probability
  eer = clamp(eer, 0, 1);

  // Calculate ARR and NNT
  const arr = cer - eer;
  const absArr = Math.abs(arr);
  const nnt = absArr < 0.0001 ? Infinity : Math.ceil(1 / absArr);
  const isHarm = arr < 0;

  // Update display
  setText('nntValue', nnt === Infinity ? '\u221E' : nnt.toString());
  setText('nntLabel', isHarm ? 'NNH (Harm)' : 'NNT (Benefit)');
  setText('arrValue', `${(arr * 100).toFixed(1)}%`);

  // Cost calculations
  const costPerEventAvoided = nnt === Infinity ? Infinity : nnt * costTreatment;
  setText('costAvoided',
    costPerEventAvoided === Infinity ? '\u221E' : `\u00A3${costPerEventAvoided.toLocaleString()}`
  );

  // ICER calculation
  const qalysGained = absArr * qalyLoss;
  const icer = qalysGained > 0 ? costTreatment / qalysGained : Infinity;

  setText('icerValue',
    icer === Infinity ? '\u221E' : `\u00A3${Math.round(icer).toLocaleString()}`
  );

  const icerThresholdEl = getElement('icerThreshold');
  if (icerThresholdEl) {
    if (icer === Infinity) {
      icerThresholdEl.textContent = 'Cannot calculate';
    } else if (icer < WTP_THRESHOLD) {
      icerThresholdEl.textContent = `\u2713 Below \u00A3${WTP_THRESHOLD.toLocaleString()} threshold`;
      icerThresholdEl.style.color = 'var(--success)';
    } else {
      icerThresholdEl.textContent = `\u2717 Above \u00A3${WTP_THRESHOLD.toLocaleString()} threshold`;
      icerThresholdEl.style.color = 'var(--warning)';
    }
  }

  // Render NNT grid
  renderNNTGrid(nnt, isHarm);

  // Update NNT info section
  const nntDisplay = getElement('nntDisplay');
  if (nntDisplay) {
    nntDisplay.textContent = nnt === Infinity ? '\u221E' : nnt.toString();
    nntDisplay.style.color = isHarm ? 'var(--warning)' : 'var(--success)';
  }

  setText('nntType', isHarm ? 'NNH' : 'NNT');

  const nntDesc = getElement('nntDesc');
  if (nntDesc) {
    if (nnt === Infinity) {
      nntDesc.textContent = 'No absolute difference between groups';
    } else {
      nntDesc.textContent = `Treat ${nnt} patients at ${(cer * 100).toFixed(0)}% baseline risk for 1 ${isHarm ? 'additional harm' : 'benefit'}`;
    }
  }
}

/**
 * Render the visual NNT grid (icon array)
 * @param {number} nnt - Number needed to treat
 * @param {boolean} isHarm - Whether effect is harmful
 */
function renderNNTGrid(nnt, isHarm) {
  const container = getElement('nntGrid');
  if (!container) return;

  const maxDots = 100;
  const activeDots = nnt === Infinity ? 0 : Math.min(nnt, maxDots);

  let html = '';
  for (let i = 0; i < maxDots; i++) {
    const isActive = i < activeDots;
    const activeClass = isActive ? (isHarm ? 'active' : 'active benefit') : '';
    html += `<div class="nnt-dot ${activeClass}" title="${i + 1}"></div>`;
  }

  container.innerHTML = html;
}

/**
 * Calculate NNT from OR and baseline risk
 * @param {number} or - Odds ratio
 * @param {number} cer - Control event rate
 * @returns {Object} NNT calculation result
 */
export function calculateNNTFromOR(or, cer) {
  if (or <= 0 || cer <= 0 || cer >= 1) {
    return { nnt: Infinity, arr: 0, eer: cer };
  }

  const eer = (cer * or) / (1 - cer + cer * or);
  const arr = cer - eer;
  const nnt = Math.abs(arr) < 0.0001 ? Infinity : Math.ceil(1 / Math.abs(arr));

  return { nnt, arr, eer, isHarm: arr < 0 };
}

/**
 * Calculate NNT from RR and baseline risk
 * @param {number} rr - Risk ratio
 * @param {number} cer - Control event rate
 * @returns {Object} NNT calculation result
 */
export function calculateNNTFromRR(rr, cer) {
  if (rr <= 0 || cer <= 0 || cer >= 1) {
    return { nnt: Infinity, arr: 0, eer: cer };
  }

  const eer = clamp(cer * rr, 0, 1);
  const arr = cer - eer;
  const nnt = Math.abs(arr) < 0.0001 ? Infinity : Math.ceil(1 / Math.abs(arr));

  return { nnt, arr, eer, isHarm: arr < 0 };
}

/**
 * Calculate ICER
 * @param {number} incrementalCost - Incremental cost of treatment
 * @param {number} incrementalQALY - Incremental QALY gained
 * @returns {number} ICER value
 */
export function calculateICER(incrementalCost, incrementalQALY) {
  if (incrementalQALY <= 0) return Infinity;
  return incrementalCost / incrementalQALY;
}
