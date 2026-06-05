# HTA Meta-Analysis Suite

A single-file, offline browser tool for health-technology-assessment style
meta-analysis. Enter study-level data (binary, continuous, single proportion,
survival HR/CI, pre-computed effect/SE, or correlation), pool on the log scale
where appropriate, and inspect results across forest, funnel, Baujat, L'Abbe,
radial, GOSH, TSA, sensitivity (leave-one-out / cumulative / influence), GRADE,
risk-of-bias, and simple health-economic views.

The statistics run entirely in the browser. No network access, no build step.

## What it computes

- **Pooling**: inverse-variance random-effects and fixed-effect, pooled on the
  log scale for ratio measures (OR, RR, HR) and back-transformed for display.
- **tau-squared estimators**: DerSimonian-Laird, REML, Paule-Mandel,
  Sidik-Jonkman (plus UI options that fall back to DL).
- **Heterogeneity**: Cochran's Q, I-squared, H-squared, tau.
- **HKSJ (Knapp-Hartung)** SE adjustment with a `max(1, q)` variance floor and a
  Student-t critical value on `k-1` degrees of freedom.
- **Prediction interval** using `t_{k-1} * sqrt(tau^2 + SE^2)`.
- **Publication bias**: Egger's regression, Begg's rank correlation,
  Rosenthal fail-safe N, trim-and-fill.
- **TSA**: required information size, diversity (D-squared), O'Brien-Fleming
  monitoring boundary, cumulative Z-curve.

## Layout

| File          | Purpose                                                              |
|---------------|---------------------------------------------------------------------|
| `index.html`  | The application (UI + DOM/Plotly rendering). Loads the two scripts.  |
| `engine.js`   | Pure statistical engine, extracted verbatim. Single source of truth. |
| `plotly.min.js` | Vendored Plotly 2.27.0 (offline; no CDN).                         |
| `tests.js`    | Pure-Node unit tests for the engine.                                |

`index.html` loads `plotly.min.js` and `engine.js` before its inline script, so
the page works fully offline (open the file directly or serve the folder).

## Running the tests

```
node tests.js
```

Prints each check and a final `N passed, M failed` line; exits non-zero on any
failure. Expected values are hand-derived independently of the engine, including
a fully worked 2-study OR pooling example and a 3-study DL example with real
heterogeneity, plus edge cases (k=1 passthrough, two identical studies giving
tau^2 = 0 / I^2 = 0, and empty / all-excluded guards).

## Fixes applied during revival

- **Offline**: vendored Plotly 2.27.0 locally as `plotly.min.js` and removed the
  Google Fonts `<link>` (replaced with a comment). Zero external URLs remain in
  `index.html` (verified: no `http(s)://`, no `@import url(...)`).
- **Engine extraction**: the pure `Engine` object was extracted verbatim into
  `engine.js` (with a Node `module.exports`) and the inline duplicate was
  deleted, so the browser and the tests share one implementation.
- **Tests added**: `tests.js` (45 checks) with hand-derived expectations.

No statistical methodology was changed. During testing the random-effects
`poolIV` was confirmed to report Q / I-squared / H-squared using random-effects
weights (which equals `(k-1)` times the Knapp-Hartung q it then uses for the
HKSJ floor), while the tau-squared estimators use the conventional
inverse-variance-weighted Cochran's Q internally; the tests document both.
