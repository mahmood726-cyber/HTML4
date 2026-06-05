// Statistical engine for HTA Meta-Analysis Suite
// Extracted VERBATIM from index.html (single source of truth). Pure functions only.

const Engine = {
    // Utility functions
    num: v => { const n = parseFloat(v); return isNaN(n) ? null : n; },
    
    // Normal CDF (Abramowitz & Stegun)
    normCDF: z => {
        const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
        const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;
        const sign = z < 0 ? -1 : 1;
        z = Math.abs(z) / Math.SQRT2;
        const t = 1 / (1 + p * z);
        const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t * Math.exp(-z*z);
        return 0.5 * (1 + sign * y);
    },
    
    // P-values
    pFromZ: z => 2 * (1 - Engine.normCDF(Math.abs(z))),
    
    // T-distribution critical value (improved approximation)
    tCrit: (df, alpha = 0.05) => {
        if (df <= 0) return 1.96;
        if (df === 1) return 12.706;
        if (df === 2) return 4.303;
        if (df === 3) return 3.182;
        if (df === 4) return 2.776;
        if (df === 5) return 2.571;
        const z = 1.96;
        return z + (z**3 + z) / (4*df) + (5*z**5 + 16*z**3 + 3*z) / (96*df**2);
    },
    
    // T-distribution CDF approximation
    tCDF: (t, df) => {
        if (df > 100) return Engine.normCDF(t);
        const z = t * Math.sqrt(1 - 1/(4*df) - 7/(120*df*df));
        return Engine.normCDF(z);
    },
    
    pFromT: (t, df) => 2 * (1 - Engine.tCDF(Math.abs(t), df)),
    
    // Chi-squared p-value
    chiSqP: (x, df) => {
        if (x <= 0 || df <= 0) return 1;
        const z = Math.pow(x/df, 1/3) - (1 - 2/(9*df));
        const se = Math.sqrt(2/(9*df));
        return 1 - Engine.normCDF(z/se);
    },

    // Effect size calculations
    calcEffect: (row, type, metric) => {
        const cc = 0.5;
        
        if (type === 'binary') {
            let e1 = Engine.num(row.e1), n1 = Engine.num(row.n1);
            let e2 = Engine.num(row.e2), n2 = Engine.num(row.n2);
            if ([e1,n1,e2,n2].some(v => v === null || v < 0)) return null;
            if (n1 === 0 || n2 === 0) return null;

            const needsCC = e1 === 0 || e2 === 0 || e1 === n1 || e2 === n2;
            const e1c = needsCC ? e1 + cc : e1;
            const n1c = needsCC ? n1 + 2*cc : n1;
            const e2c = needsCC ? e2 + cc : e2;
            const n2c = needsCC ? n2 + 2*cc : n2;

            const p1 = e1c/n1c, p2 = e2c/n2c;
            let es, vi, display;

            if (metric === 'OR') {
                es = Math.log((e1c*(n2c-e2c)) / (e2c*(n1c-e1c)));
                vi = 1/e1c + 1/(n1c-e1c) + 1/e2c + 1/(n2c-e2c);
                display = Math.exp(es);
            } else if (metric === 'RR') {
                es = Math.log(p1/p2);
                vi = (1-p1)/e1c + (1-p2)/e2c;
                display = Math.exp(es);
            } else { // RD
                es = p1 - p2;
                vi = p1*(1-p1)/n1c + p2*(1-p2)/n2c;
                display = es;
            }

            return { es, vi, se: Math.sqrt(vi), display, raw: {e1,n1,e2,n2,p1:e1/n1,p2:e2/n2} };
        }

        if (type === 'continuous') {
            const m1 = Engine.num(row.m1), s1 = Engine.num(row.s1), n1 = Engine.num(row.n1);
            const m2 = Engine.num(row.m2), s2 = Engine.num(row.s2), n2 = Engine.num(row.n2);
            if ([m1,s1,n1,m2,s2,n2].some(v => v === null)) return null;
            if (s1 <= 0 || s2 <= 0 || n1 <= 0 || n2 <= 0) return null;

            const pooledSD = Math.sqrt(((n1-1)*s1*s1 + (n2-1)*s2*s2) / (n1+n2-2));
            const d = (m1 - m2) / pooledSD;
            const j = 1 - 3/(4*(n1+n2-2) - 1);
            const es = d * j;
            const vi = (n1+n2)/(n1*n2) + (es*es)/(2*(n1+n2));

            return { es, vi, se: Math.sqrt(vi), display: es };
        }

        if (type === 'proportion') {
            const e = Engine.num(row.e), n = Engine.num(row.n);
            if (e === null || n === null || n <= 0) return null;

            const p_adj = (e + 0.5) / (n + 1);
            const es = Math.log(p_adj / (1 - p_adj));
            const vi = 1 / (n * p_adj * (1 - p_adj));

            return { es, vi, se: Math.sqrt(vi), display: e/n, raw: {e, n} };
        }

        if (type === 'survival') {
            const hr = Engine.num(row.hr);
            const ll = Engine.num(row.ll), ul = Engine.num(row.ul);
            if (hr === null || hr <= 0) return null;

            const es = Math.log(hr);
            let vi;
            if (ll && ul && ll > 0 && ul > 0) {
                const se = (Math.log(ul) - Math.log(ll)) / (2 * 1.96);
                vi = se * se;
            } else {
                const se = Engine.num(row.se);
                if (!se || se <= 0) return null;
                vi = se * se;
            }

            return { es, vi, se: Math.sqrt(vi), display: hr };
        }

        if (type === 'correlation') {
            const r = Engine.num(row.r), n = Engine.num(row.n);
            if (r === null || n === null || n <= 3) return null;
            if (Math.abs(r) >= 1) return null;

            // Fisher's z transformation
            const es = 0.5 * Math.log((1 + r) / (1 - r));
            const vi = 1 / (n - 3);

            return { es, vi, se: Math.sqrt(vi), display: r };
        }

        if (type === 'generic') {
            const es = Engine.num(row.es), se = Engine.num(row.se);
            if (es === null || se === null || se <= 0) return null;
            return { es, vi: se*se, se, display: es };
        }

        return null;
    },

    // Tau² estimators
    tauDL: (effects) => {
        const k = effects.length;
        if (k < 2) return 0;
        const w = effects.map(e => 1/e.vi);
        const sumW = w.reduce((a,b) => a+b, 0);
        const sumW2 = w.reduce((a,b) => a + b*b, 0);
        const mu = effects.reduce((a,e,i) => a + w[i]*e.es, 0) / sumW;
        const Q = effects.reduce((a,e,i) => a + w[i] * (e.es - mu)**2, 0);
        const C = sumW - sumW2/sumW;
        return Math.max(0, (Q - (k-1)) / C);
    },

    tauREML: (effects, maxIter = 50, tol = 1e-6) => {
        let tau2 = Engine.tauDL(effects);
        const k = effects.length;
        for (let i = 0; i < maxIter; i++) {
            const w = effects.map(e => 1/(e.vi + tau2));
            const sumW = w.reduce((a,b) => a+b, 0);
            const mu = effects.reduce((a,e,j) => a + w[j]*e.es, 0) / sumW;
            const num = effects.reduce((a,e,j) => a + w[j]**2 * ((e.es-mu)**2 - e.vi), 0);
            const den = effects.reduce((a,e,j) => a + w[j]**2, 0);
            const tau2_new = Math.max(0, tau2 + num/den);
            if (Math.abs(tau2_new - tau2) < tol) break;
            tau2 = tau2_new;
        }
        return tau2;
    },

    tauPM: (effects, maxIter = 100, tol = 1e-6) => {
        let tau2 = Engine.tauDL(effects);
        const k = effects.length;
        for (let i = 0; i < maxIter; i++) {
            const w = effects.map(e => 1/(e.vi + tau2));
            const sumW = w.reduce((a,b) => a+b, 0);
            const mu = effects.reduce((a,e,j) => a + w[j]*e.es, 0) / sumW;
            const Q = effects.reduce((a,e,j) => a + w[j] * (e.es - mu)**2, 0);
            if (Math.abs(Q - (k-1)) < tol) break;
            const dQ = -effects.reduce((a,e,j) => a + w[j]**2 * ((e.es - mu)**2 - (1 - w[j]/sumW)/w[j]), 0);
            if (Math.abs(dQ) < tol) break;
            tau2 = Math.max(0, tau2 - (Q - (k-1))/dQ);
        }
        return tau2;
    },

    tauSJ: (effects) => {
        const k = effects.length;
        if (k < 2) return 0;
        const mu0 = effects.reduce((a,e) => a + e.es, 0) / k;
        const tau2_0 = effects.reduce((a,e) => a + (e.es - mu0)**2, 0) / (k-1);
        const w = effects.map(e => 1/(e.vi + tau2_0));
        const sumW = w.reduce((a,b) => a+b, 0);
        const mu = effects.reduce((a,e,j) => a + w[j]*e.es, 0) / sumW;
        const Q = effects.reduce((a,e,j) => a + w[j] * (e.es - mu)**2, 0);
        return Math.max(0, (Q - (k-1)) * tau2_0 / Q);
    },

    // Pool effects
    poolIV: (effects, tau2, useHKSJ = true, confLevel = 0.95) => {
        const active = effects.filter(e => !e.excluded);
        if (active.length === 0) return null;

        const w = active.map(e => 1/(e.vi + tau2));
        const sumW = w.reduce((a,b) => a+b, 0);
        const es = active.reduce((a,e,i) => a + w[i]*e.es, 0) / sumW;
        let se = Math.sqrt(1/sumW);

        const Q = active.reduce((a,e,i) => a + w[i] * (e.es - es)**2, 0);
        const df = active.length - 1;
        const I2 = df > 0 ? Math.max(0, (Q - df) / Q * 100) : 0;
        const H2 = df > 0 ? Q / df : 1;

        let tCrit = 1.96;
        const alpha = 1 - confLevel;
        if (useHKSJ && df > 0) {
            const qStar = Q / df;
            se = se * Math.sqrt(Math.max(1, qStar));
            tCrit = Engine.tCrit(df, alpha);
        }

        const z = es / se;
        const pVal = useHKSJ && df > 0 ? Engine.pFromT(z, df) : Engine.pFromZ(z);
        const ciLo = es - tCrit * se;
        const ciHi = es + tCrit * se;

        // Prediction interval
        const piSe = Math.sqrt(se*se + tau2);
        const piCrit = df > 0 ? Engine.tCrit(df, alpha) : 1.96;
        const piLo = es - piCrit * piSe;
        const piHi = es + piCrit * piSe;

        // Study weights and contributions
        effects.forEach(e => {
            if (e.excluded) {
                e.w = 0; e.q = 0; e.imp = 0;
            } else {
                const wi = 1/(e.vi + tau2);
                e.w = wi / sumW * 100;
                e.q = wi * (e.es - es)**2;
                const sumW_loo = sumW - wi;
                const es_loo = (es * sumW - wi * e.es) / sumW_loo;
                e.imp = Math.abs(es - es_loo);
            }
        });

        return {
            es, se, z, pVal,
            ciLo, ciHi,
            piLo, piHi,
            Q, df, I2, H2, tau2,
            k: active.length,
            studies: effects,
            totalN: active.reduce((a,s) => {
                if (s.raw) return a + (s.raw.n1 || 0) + (s.raw.n2 || 0) + (s.raw.n || 0);
                return a;
            }, 0)
        };
    },

    // Publication bias tests
    eggerTest: (effects) => {
        const active = effects.filter(e => !e.excluded);
        const n = active.length;
        if (n < 3) return { p: null, intercept: null };

        const x = active.map(e => 1/e.se);
        const y = active.map(e => e.es/e.se);
        const sumX = x.reduce((a,b) => a+b, 0);
        const sumY = y.reduce((a,b) => a+b, 0);
        const sumXY = x.reduce((a,xi,i) => a + xi*y[i], 0);
        const sumX2 = x.reduce((a,xi) => a + xi*xi, 0);

        const slope = (n*sumXY - sumX*sumY) / (n*sumX2 - sumX*sumX);
        const intercept = (sumY - slope*sumX) / n;

        const yPred = x.map(xi => intercept + slope*xi);
        const sse = y.reduce((a,yi,i) => a + (yi - yPred[i])**2, 0);
        const mse = sse / (n-2);
        const sxx = sumX2 - sumX*sumX/n;
        const seInt = Math.sqrt(mse * (1/n + (sumX/n)**2/sxx));

        const t = intercept / seInt;
        const p = Engine.pFromT(t, n-2);

        return { p, intercept, t };
    },

    beggTest: (effects) => {
        const active = effects.filter(e => !e.excluded);
        const n = active.length;
        if (n < 3) return { p: null, tau: null };

        // Rank correlation (Kendall's tau)
        const sorted = [...active].sort((a,b) => a.es - b.es);
        const ranks = sorted.map((s,i) => ({...s, rank: i+1}));
        
        let concordant = 0, discordant = 0;
        for (let i = 0; i < n-1; i++) {
            for (let j = i+1; j < n; j++) {
                const seComp = ranks[i].se - ranks[j].se;
                const rankComp = ranks[i].rank - ranks[j].rank;
                if (seComp * rankComp > 0) concordant++;
                else if (seComp * rankComp < 0) discordant++;
            }
        }
        
        const tau = (concordant - discordant) / (n * (n-1) / 2);
        const z = 3 * tau * Math.sqrt(n * (n-1)) / Math.sqrt(2 * (2*n + 5));
        const p = Engine.pFromZ(z);

        return { p, tau };
    },

    failSafeN: (effects) => {
        const active = effects.filter(e => !e.excluded);
        const zCrit = 1.645;
        const zSum = active.reduce((a,e) => a + e.es/e.se, 0);
        return Math.max(0, Math.floor((zSum*zSum)/(zCrit*zCrit) - active.length));
    },

    trimAndFill: (effects) => {
        const active = effects.filter(e => !e.excluded).map(e => ({...e}));
        const n = active.length;
        if (n < 3) return { k0: 0, adjusted: null };

        active.sort((a,b) => a.es - b.es);
        const median = active[Math.floor(n/2)].es;
        const left = active.filter(e => e.es < median).length;
        const right = active.filter(e => e.es > median).length;
        const k0 = Math.abs(right - left);

        if (k0 > 0) {
            const toMirror = right > left ? active.slice(-k0) : active.slice(0, k0);
            const imputed = toMirror.map(e => ({
                es: 2*median - e.es,
                vi: e.vi,
                se: e.se,
                excluded: false,
                imputed: true
            }));
            const augmented = [...active, ...imputed];
            const tau2 = Engine.tauDL(augmented);
            const adjusted = Engine.poolIV(augmented, tau2, false);
            return { k0, adjusted };
        }
        return { k0: 0, adjusted: null };
    },

    // TSA calculations
    calcTSA: (effects, alpha = 0.05, beta = 0.20, rrr = 0.20, cer = 0.10) => {
        const active = effects.filter(e => !e.excluded);
        if (active.length < 2) return null;

        // Calculate diversity (D²) - adjusted heterogeneity
        const tau2 = Engine.tauDL(active);
        const avgVi = active.reduce((a,e) => a + e.vi, 0) / active.length;
        const D2 = tau2 / (tau2 + avgVi) * 100;

        // Required information size (simplified)
        const zAlpha = 1.96; // two-sided
        const zBeta = 0.84; // 80% power
        const pC = cer;
        const pE = cer * (1 - rrr);
        const pooledP = (pC + pE) / 2;
        
        // Sample size per group for single trial
        const nPerGroup = 2 * pooledP * (1 - pooledP) * (zAlpha + zBeta)**2 / (pC - pE)**2;
        
        // Diversity-adjusted RIS
        const RIS = Math.ceil(2 * nPerGroup / (1 - D2/100));

        // Accrued information
        const accrued = active.reduce((a,s) => {
            if (s.raw) return a + (s.raw.n1 || 0) + (s.raw.n2 || 0);
            return a;
        }, 0);

        // Cumulative Z-values for TSA plot
        const cumZ = [];
        const cumN = [];
        let cumES = 0, cumW = 0;
        
        active.forEach((s, i) => {
            const w = 1 / s.vi;
            cumW += w;
            cumES = (cumES * (cumW - w) + w * s.es) / cumW;
            const cumSE = Math.sqrt(1 / cumW);
            cumZ.push(cumES / cumSE);
            cumN.push(active.slice(0, i+1).reduce((a,x) => {
                if (x.raw) return a + (x.raw.n1 || 0) + (x.raw.n2 || 0);
                return a + 100; // default if no raw data
            }, 0));
        });

        // O'Brien-Fleming spending function for monitoring boundary
        const infoFracs = cumN.map(n => n / RIS);
        const monitoringBoundary = infoFracs.map(t => {
            if (t <= 0) return 8;
            return zAlpha / Math.sqrt(t);
        });

        // Futility boundary (simplified)
        const futilityBoundary = infoFracs.map(t => {
            if (t <= 0) return 0;
            return -0.5 + zBeta * Math.sqrt(t);
        });

        // Determine conclusion
        let conclusion = 'Inconclusive';
        let boundary = 'Continue accrual';
        const lastZ = cumZ[cumZ.length - 1];
        const lastMonitor = monitoringBoundary[monitoringBoundary.length - 1];
        
        if (Math.abs(lastZ) >= lastMonitor) {
            conclusion = lastZ > 0 ? 'Benefit' : 'Harm';
            boundary = 'Crossed monitoring boundary';
        } else if (accrued >= RIS * 0.9) {
            conclusion = 'Futility likely';
            boundary = 'Near RIS without crossing';
        }

        return {
            RIS,
            accrued,
            pctRIS: (accrued / RIS * 100).toFixed(1),
            D2: D2.toFixed(1),
            cumZ,
            cumN,
            monitoringBoundary,
            futilityBoundary,
            conclusion,
            boundary
        };
    }
};

if (typeof module!=="undefined"&&module.exports){ module.exports = Engine; }
