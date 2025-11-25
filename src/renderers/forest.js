/**
 * Forest plot renderer
 */

import { escapeHTML, formatEffectSize } from '../utils/index.js';
import { getElement } from '../app/dom.js';

/**
 * Render the forest plot
 * @param {Object} result - Pooled meta-analysis result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderForest(result, isRatio) {
  const container = getElement('forestPlot');
  if (!container) return;

  const active = result.studies.filter(s => !s.excluded);
  if (active.length === 0) {
    container.innerHTML = '<p class="plot-message">No studies to display</p>';
    return;
  }

  // Calculate plot bounds
  const allVals = active.flatMap(s => [s.es - 1.96 * s.se, s.es + 1.96 * s.se]);
  allVals.push(result.ciLo, result.ciHi);

  let min = Math.min(...allVals);
  let max = Math.max(...allVals);
  const pad = (max - min) * 0.15;
  min -= pad;
  max += pad;

  const toPos = v => ((v - min) / (max - min) * 100);
  const fmtVal = v => formatEffectSize(v, isRatio);

  // Build HTML
  let html = `
    <div class="forest-row header">
      <div class="forest-study">Study</div>
      <div class="forest-weight">Weight</div>
      <div class="forest-plot-area"></div>
      <div class="forest-es">Effect</div>
      <div class="forest-ci">95% CI</div>
    </div>
  `;

  // Study rows
  active.forEach(s => {
    const lo = s.es - 1.96 * s.se;
    const hi = s.es + 1.96 * s.se;
    const posLo = toPos(lo);
    const posHi = toPos(hi);
    const posES = toPos(s.es);

    // Scale point size by weight
    const size = Math.max(6, Math.min(14, 6 + s.w * 0.25));

    const yearStr = s.year ? ` (${escapeHTML(s.year)})` : '';

    html += `
      <div class="forest-row">
        <div class="forest-study">${escapeHTML(s.id)}${yearStr}</div>
        <div class="forest-weight">${s.w.toFixed(1)}%</div>
        <div class="forest-plot-area">
          <div class="forest-null-line" style="left:${toPos(0)}%"></div>
          <div class="forest-ci-line" style="left:${posLo}%;width:${posHi - posLo}%"></div>
          <div class="forest-point" style="left:${posES}%;width:${size}px;height:${size}px"></div>
        </div>
        <div class="forest-es">${fmtVal(s.es)}</div>
        <div class="forest-ci">[${fmtVal(lo)}, ${fmtVal(hi)}]</div>
      </div>
    `;
  });

  // Summary row
  const sLo = toPos(result.ciLo);
  const sHi = toPos(result.ciHi);
  const sES = toPos(result.es);

  html += `
    <div class="forest-row summary">
      <div class="forest-study">RE Model (k=${result.k})</div>
      <div class="forest-weight">100%</div>
      <div class="forest-plot-area">
        <div class="forest-null-line" style="left:${toPos(0)}%;background:rgba(255,255,255,0.3)"></div>
        <div class="forest-ci-line" style="left:${sLo}%;width:${sHi - sLo}%;background:var(--paper)"></div>
        <div class="forest-point" style="left:${sES}%"></div>
      </div>
      <div class="forest-es">${fmtVal(result.es)}</div>
      <div class="forest-ci">[${fmtVal(result.ciLo)}, ${fmtVal(result.ciHi)}]</div>
    </div>
  `;

  container.innerHTML = html;
}

/**
 * Render leave-one-out forest plot
 * @param {Object} result - Pooled result
 * @param {boolean} isRatio - Whether effect is ratio-based
 */
export function renderLOO(result, isRatio) {
  const container = getElement('looPlot');
  if (!container) return;

  const active = result.studies.filter(s => !s.excluded);
  if (active.length < 3) {
    container.innerHTML = '<p class="plot-message">Need ≥3 studies for leave-one-out analysis</p>';
    return;
  }

  // Import dynamically to avoid circular dependencies
  import('../engine/tau.js').then(({ tauDL }) => {
    import('../engine/pooling.js').then(({ poolInverseVariance }) => {
      const looResults = [];

      active.forEach((excluded, i) => {
        const subset = active.filter((_, j) => j !== i);
        const tau2 = tauDL(subset);
        const pooled = poolInverseVariance(
          subset.map(s => ({ ...s })),
          tau2,
          { useHKSJ: false }
        );

        if (pooled) {
          looResults.push({
            label: `Excl. ${excluded.id}`,
            es: pooled.es,
            lo: pooled.ciLo,
            hi: pooled.ciHi
          });
        }
      });

      // Calculate bounds
      const allVals = looResults.flatMap(r => [r.lo, r.hi]);
      allVals.push(result.ciLo, result.ciHi);

      let min = Math.min(...allVals);
      let max = Math.max(...allVals);
      const pad = (max - min) * 0.15;
      min -= pad;
      max += pad;

      const toPos = v => ((v - min) / (max - min) * 100);
      const fmtVal = v => formatEffectSize(v, isRatio);

      let html = `
        <div class="forest-row header">
          <div class="forest-study">Omitting</div>
          <div class="forest-weight"></div>
          <div class="forest-plot-area"></div>
          <div class="forest-es">ES</div>
          <div class="forest-ci">95% CI</div>
        </div>
      `;

      looResults.forEach(r => {
        html += `
          <div class="forest-row">
            <div class="forest-study">${escapeHTML(r.label)}</div>
            <div class="forest-weight"></div>
            <div class="forest-plot-area">
              <div class="forest-null-line" style="left:${toPos(0)}%"></div>
              <div class="forest-ci-line" style="left:${toPos(r.lo)}%;width:${toPos(r.hi) - toPos(r.lo)}%"></div>
              <div class="forest-point" style="left:${toPos(r.es)}%"></div>
            </div>
            <div class="forest-es">${fmtVal(r.es)}</div>
            <div class="forest-ci">[${fmtVal(r.lo)}, ${fmtVal(r.hi)}]</div>
          </div>
        `;
      });

      // Overall summary
      html += `
        <div class="forest-row summary">
          <div class="forest-study">All Studies</div>
          <div class="forest-weight"></div>
          <div class="forest-plot-area">
            <div class="forest-null-line" style="left:${toPos(0)}%;background:rgba(255,255,255,0.3)"></div>
            <div class="forest-ci-line" style="left:${toPos(result.ciLo)}%;width:${toPos(result.ciHi) - toPos(result.ciLo)}%;background:var(--paper)"></div>
            <div class="forest-point" style="left:${toPos(result.es)}%"></div>
          </div>
          <div class="forest-es">${fmtVal(result.es)}</div>
          <div class="forest-ci">[${fmtVal(result.ciLo)}, ${fmtVal(result.ciHi)}]</div>
        </div>
      `;

      container.innerHTML = html;
    });
  });
}
