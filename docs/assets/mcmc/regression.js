/* Small, genuine von Mises regression models; existing examples stay unchanged. */
(function(root) {
  'use strict';
  const C = typeof module === 'object' && module.exports ? require('./core.js') : root
    .CircularInference;

  function validate(data, design, hierarchical) {
    if (!data.length || design.x?.length !== data.length || data.some(x => !Number.isFinite(x)) ||
      design.x.some(x => !Number.isFinite(x))) throw Error('Invalid regression data');
    if (hierarchical && (design.individual?.length !== data.length || design.individual.some(i => !
        Number.isInteger(i) || i < 0 || i >= 4))) throw Error('Invalid individual indices');
  }

  function regression(data, design) {
    validate(data, design, false);
    const x = design.x,
      n = data.length,
      vary = !!design.varyConcentration;

    function target(q) {
      const [a, b, e, g = 0] = q, grad = [-a / C.TAU ** 2, -b / 2.25, -e / 9, ...(vary ? [-g / 2.25] :
        [])];
      let logp = -.5 * (a / C.TAU) ** 2 - .5 * (b / 1.5) ** 2 - .5 * (e / 3) ** 2 - (vary ? .5 * (g /
        1.5) ** 2 : 0);
      data.forEach((y, i) => {
        const eta = e + (vary ? g * x[i] : 0),
          k = C.softplus(eta),
          B = C.bessel(k),
          delta = y - a - b * x[i],
          gm = k * Math.sin(delta),
          ge = C.sigmoid(eta) * (Math.cos(delta) - B.A);
        logp += k * Math.cos(delta) - B.logI0;
        grad[0] += gm;
        grad[1] += x[i] * gm;
        grad[2] += ge;
        if (vary) grad[3] += x[i] * ge;
      });
      return {
        logp,
        grad
      };
    }
    let best = null;
    for (let j = -60; j <= 60; j++) {
      const slope = j * Math.PI / 60,
        fit = C.model(data.map((y, i) => C.wrap(y - slope * x[i])));
      if (!best || fit.R > best.fit.R) best = {
        slope,
        fit
      };
    }
    return {
      kind: 'regression',
      n,
      khat: best.fit.khat,
      target,
      components(q, at = {}) {
        const x = at.x || 0,
          eta = q[0] + q[1] * x,
          etaK = q[2] + (vary ? q[3] * x : 0);
        return [{
          index: 0,
          mu: C.wrap(eta),
          etaMu: eta,
          etaK,
          k: C.softplus(etaK),
          weight: 1
        }];
      },
      initialize(random, index) {
        const k = Math.min(12, best.fit.khat),
          e = C.invsoftplus(k),
          information = Math.max(1, n * k * best.fit.R);
        return {
          q: [best.fit.mean + [-C.TAU, 0, C.TAU, 0][index % 4] + .1 * C.normal(random), best.slope +
            .1 * C.normal(random), e + .2 * C.normal(random), ...(vary ? [.1 * C.normal(random)] :
            [])
          ],
          mass: [information, Math.max(1, information * x.reduce((s, v) => s + v * v, 0) / n), Math
            .max(1, n * .1), ...(vary ? [Math.max(1, n * .04)] : [])
          ]
        };
      }
    };
  }

  function hierarchical(data, design) {
    validate(data, design, true);
    const n = data.length,
      x = design.x,
      ids = design.individual,
      parts = [0, 1].map(t => C.model(data.filter((_, i) => x[i] === t)));

    function target(q) {
      const [a, b, e, l] = q, tau = C.softplus(l), dt = C.sigmoid(l), sigma = C.softplus(q[9]), ds = C
        .sigmoid(q[9]);
      let logp = -.5 * (a / C.TAU) ** 2 - .5 * (b / 1.5) ** 2 - .5 * (e / 3) ** 2 - .5 * (q[8] / 3) **
        2 - .5 * (tau / .8) ** 2 - Math.log1p(Math.exp(-l)) - .5 * sigma * sigma - Math.log1p(Math.exp(
          -q[9]));
      const grad = [-a / C.TAU ** 2, -b / 2.25, -e / 9, -tau * dt / .64 + 1 - dt, ...q.slice(4, 8).map(
        z => -z), -q[8] / 9, -sigma * ds + 1 - ds, ...q.slice(10, 14).map(z => -z)];
      for (const z of [...q.slice(4, 8), ...q.slice(10, 14)]) logp -= .5 * z * z;
      data.forEach((y, i) => {
        const t = x[i],
          id = ids[i],
          ki = t ? 8 : 2,
          v = q[10 + id],
          etaK = q[ki] + sigma * v,
          k = C.softplus(etaK),
          B = C.bessel(k),
          z = q[4 + id],
          delta = y - a - b * t - tau * z,
          gm = k * Math.sin(delta),
          gk = C.sigmoid(etaK) * (Math.cos(delta) - B.A);
        logp += k * Math.cos(delta) - B.logI0;
        grad[0] += gm;
        grad[1] += t * gm;
        grad[ki] += gk;
        grad[3] += gm * dt * z;
        grad[4 + id] += gm * tau;
        grad[9] += gk * ds * v;
        grad[10 + id] += gk * sigma;
      });
      return {
        logp,
        grad
      };
    }
    return {
      kind: 'hierarchical',
      targetAcceptance: .95,
      n,
      khat: Math.max(...parts.map(p => p.khat)),
      target,
      components(q, at = {}) {
        const offset = at.individual >= 0 ? C.softplus(q[3]) * q[4 + at.individual] : 0,
          kOffset = at.individual >= 0 ? C.softplus(q[9]) * q[10 + at.individual] : 0;
        return [0, 1].map(t => {
          const eta = q[0] + q[1] * t + offset,
            etaK = q[t ? 8 : 2] + kOffset;
          return {
            index: t,
            mu: C.wrap(eta),
            etaMu: eta,
            etaK,
            k: C.softplus(etaK),
            weight: 1
          };
        });
      },
      initialize(random, index) {
        const a = parts[0].mean,
          b = C.wrap(parts[1].mean - a),
          ks = parts.map(p => Math.min(10, Math.max(1, p.khat))),
          k = (ks[0] + ks[1]) / 2,
          tau = .4,
          information = Math.max(1, n * k * .7);
        const z = Array.from({
          length: 4
        }, (_, j) => {
          const fit = C.model(data.filter((_, i) => ids[i] === j));
          return C.wrap(fit.mean - a - b * .5) / tau;
        });
        return {
          q: [a + [-C.TAU, 0, C.TAU, 0][index % 4] + .08 * C.normal(random), b + .08 * C.normal(
            random), C.invsoftplus(ks[0]) + .15 * C.normal(random), C.invsoftplus(tau) + .1 * C
            .normal(random), ...z.map(v => v + .1 * C.normal(random)), C.invsoftplus(ks[1]) + .15 * C
            .normal(random), C.invsoftplus(.6) + .1 * C.normal(random), ...Array.from({
              length: 4
            }, () => .05 * C.normal(random))
          ],
          mass: [information, information * .5, Math.max(1, n * .1), 16, ...Array(4).fill(Math.max(1,
            information * tau * tau / 4)), Math.max(1, n * .1), 16, ...Array(4).fill(2)]
        };
      }
    };
  }

  function generatePredictorData(kind, {
    n = 80,
    mu = 0,
    slope = 1,
    k = 4,
    k2 = k,
    spread = .35,
    kSpread = .8,
    varyConcentration = false,
    concentrationSlope = 0,
    seed = 42
  } = {}) {
    const random = C.rng(seed),
      x = [],
      individual = [],
      groups = [],
      data = [];
    let offsets = Array.from({
      length: 4
    }, () => C.normal(random) * spread);
    const average = offsets.reduce((a, b) => a + b, 0) / 4;
    offsets = offsets.map(v => v - average);
    let kOffsets = kind === 'hierarchical' ? Array.from({
      length: 4
    }, () => C.normal(random) * kSpread) : [];
    if (kOffsets.length) {
      const mean = kOffsets.reduce((a, b) => a + b, 0) / 4;
      kOffsets = kOffsets.map(v => v - mean);
    }
    for (let i = 0; i < n; i++) {
      const id = Math.floor(i / 2) % 4,
        predictor = kind === 'regression' ? -1 + 2 * i / Math.max(1, n - 1) : i % 2;
      x.push(predictor);
      individual.push(id);
      groups.push(kind === 'hierarchical' ? predictor : 0);
      data.push(C.vonMises(mu + slope * predictor + (kind === 'hierarchical' ? offsets[id] : 0),
        kind === 'hierarchical' ? C.softplus(C.invsoftplus(predictor ? k2 : k) + kOffsets[id]) :
        kind === 'regression' && varyConcentration ? C.softplus(C.invsoftplus(k) +
          concentrationSlope * predictor) : k, random));
    }
    return {
      data,
      groups,
      design: {
        x,
        individual,
        ...(kind === 'regression' ? {
          varyConcentration
        } : {})
      },
      offsets,
      kOffsets
    };
  }
  const original = C.example;
  C.example = function(data, kind, groups, design) {
    return kind === 'regression' ? regression(data, design) : kind === 'hierarchical' ? hierarchical(
      data, design) : original(data, kind, groups);
  };
  C.regression = regression;
  C.hierarchical = hierarchical;
  C.generatePredictorData = generatePredictorData;
  if (typeof module === 'object' && module.exports) module.exports = C;
})(typeof self !== 'undefined' ? self : globalThis);
