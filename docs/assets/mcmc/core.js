/* Circular likelihoods for one distribution, a mixture and two treatments.
   Priors are defined on the predictors;
   transformed densities include the pushforward Jacobian and winding sum. */
(function(root) {
  'use strict';
  const TAU = 2 * Math.PI,
    wrap = x => ((x + Math.PI) % TAU + TAU) % TAU - Math.PI;
  const softplus = x => Math.max(0, x) + Math.log1p(Math.exp(-Math.abs(x)));
  const sigmoid = x => x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x));
  const invsoftplus = x => x > 30 ? x : Math.log(Math.expm1(x));

  function rng(seed) {
    let a = seed >>> 0;
    return () => {
      a += 0x6D2B79F5;
      let t = a;
      t = Math.imul(t ^ t >>> 15, t | 1);
      t ^= t + Math.imul(t ^ t >>> 7, t | 61);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }

  function normal(random) {
    return Math.sqrt(-2 * Math.log(Math.max(1e-15, random()))) * Math.cos(TAU * random());
  }
  // Positive-term Bessel series, with its derivative, avoids cancellation for
  // the moderate concentrations in this demo. Asymptotics protect HMC tails.
  function bessel(k) {
    if (k > 50) {
      const t = 1 / k;
      return {
        logI0: k - .5 * Math.log(TAU * k) + Math.log(1 + t / 8 + 9 * t * t / 128 + 225 * t * t * t /
          3072),
        A: 1 - .5 * t - .125 * t * t - .125 * t * t * t - 25 * t * t * t * t / 128
      };
    }
    let term = 1,
      sum = 1,
      weighted = 0;
    for (let j = 1; j < 180; j++) {
      term *= k * k / (4 * j * j);
      sum += term;
      weighted += 2 * j * term;
      if (term < sum * 1e-16) break;
    }
    return {
      logI0: Math.log(sum),
      A: k > 1e-12 ? weighted / (k * sum) : k / 2
    };
  }

  function vonMises(mu, k, random) {
    if (k < 1e-6) return wrap(mu + TAU * (random() - .5));
    const a = 1 + Math.sqrt(1 + 4 * k * k),
      b = 2 * k / (a + Math.sqrt(2 * a)),
      r = (1 + b * b) / (2 * b);
    for (let i = 0; i < 10000; i++) {
      const z = Math.cos(Math.PI * random()),
        f = (1 + r * z) / (r + z),
        c = k * (r - f),
        u = Math.max(random(), 1e-15);
      if (c * (2 - c) > u || Math.log(c / u) + 1 - c >= 0) return wrap(mu + (random() > .5 ? 1 : -1) *
        Math.acos(Math.max(-1, Math.min(1, f))));
    }
    throw Error('Data generator did not converge');
  }

  function model(data) {
    const C = data.reduce((s, x) => s + Math.cos(x), 0),
      S = data.reduce((s, x) => s + Math.sin(x), 0),
      n = data.length;

    function target(q) {
      const [m, e] = q, k = softplus(e), b = bessel(k), projection = C * Math.cos(m) + S * Math.sin(m);
      return {
        logp: k * projection - n * b.logI0 - .5 * (m / TAU) ** 2 - .5 * (e / 3) ** 2,
        grad: [k * (-C * Math.sin(m) + S * Math.cos(m)) - m / (TAU * TAU), sigmoid(e) * (projection -
          n * b.A) - e / 9]
      };
    }

    function constrained(mu, k) {
      if (k <= 0) return -Infinity;
      const e = invsoftplus(k);
      let winding = 0;
      for (let j = -7; j <= 7; j++) winding += Math.exp(-.5 * ((mu + j * TAU) / TAU) ** 2);
      return k * (C * Math.cos(mu) + S * Math.sin(mu)) - n * bessel(k).logI0 - .5 * (e / 3) ** 2 + Math
        .log(winding) - Math.log(sigmoid(e));
    }
    const R = Math.hypot(C, S) / n;
    const estimate = R < .53 ? 2 * R + R ** 3 + 5 * R ** 5 / 6 : R < .85 ? -.4 + 1.39 * R + .43 / (1 -
      R) : 1 / (R ** 3 - 4 * R * R + 3 * R);
    return {
      target,
      constrained,
      n,
      C,
      S,
      R,
      mean: Math.atan2(S, C),
      khat: Math.min(30, Math.max(.1, estimate))
    };
  }

  function leapfrog(q, p, epsilon, steps, mass, target) {
    q = q.slice();
    p = p.slice();
    let t = target(q);
    const path = [q.slice()];
    for (let j = 0; j < q.length; j++) p[j] += .5 * epsilon * t.grad[j];
    for (let i = 0; i < steps; i++) {
      for (let j = 0; j < q.length; j++) q[j] += epsilon * p[j] / mass[j];
      t = target(q);
      path.push(q.slice());
      for (let j = 0; j < q.length; j++) p[j] += (i === steps - 1 ? .5 : 1) * epsilon * t.grad[j];
      if (!Number.isFinite(t.logp) || t.grad.some(x => !Number.isFinite(x))) break;
    }
    return {
      q,
      p,
      path,
      logp: t.logp
    };
  }

  function chain(m, seed, index = 0) {
    if (m.initialize) {
      const random = rng(seed),
        initial = m.initialize(random, index);
      return {
        ...initial,
        random,
        epsilon: .35,
        iteration: 0,
        accepted: 0,
        retained: 0,
        divergences: 0,
        hbar: 0,
        logAverage: Math.log(.35),
        center: Math.log(3.5)
      };
    }
    const random = rng(seed),
      k = m.khat,
      e = invsoftplus(k),
      A = bessel(k).A;
    const mass = [Math.max(1, m.n * k * m.R), Math.max(.2, m.n * Math.max(1e-4, 1 - A * A - A / k) *
      sigmoid(e) ** 2 + 1 / 9)];
    const q = [m.mean + [-TAU, 0, TAU, 0][index % 4] + .4 * normal(random), e + .5 * normal(random)];
    return {
      q,
      mass,
      random,
      epsilon: .35,
      iteration: 0,
      accepted: 0,
      retained: 0,
      divergences: 0,
      hbar: 0,
      logAverage: Math.log(.35),
      center: Math.log(3.5)
    };
  }

  function step(c, m, warmup = 120) {
    const p = c.mass.map(v => Math.sqrt(v) * normal(c.random)),
      eps = c.epsilon * (.85 + .3 * c.random()),
      L = 8 + Math.floor(17 * c.random());
    const proposal = leapfrog(c.q, p, eps, L, c.mass, m.target),
      start = m.target(c.q).logp;
    const kinetic = p => p.reduce((s, v, j) => s + v * v / (2 * c.mass[j]), 0);
    const delta = start - proposal.logp + kinetic(proposal.p) - kinetic(p);
    const alpha = Number.isFinite(delta) ? Math.exp(Math.min(0, -delta)) : 0,
      accepted = c.random() < alpha;
    const warming = c.iteration < warmup;
    if (accepted) c.q = proposal.q;
    c.iteration++;
    if (warming) {
      const t = c.iteration,
        w = 1 / (t + 10);
      c.hbar = (1 - w) * c.hbar + w * ((m.targetAcceptance || .8) - alpha);
      const logE = c.center - Math.sqrt(t) / .05 * c.hbar,
        weight = t ** (-.75);
      c.logAverage = weight * logE + (1 - weight) * c.logAverage;
      c.epsilon = Math.max(.002, Math.min(2, Math.exp(t === warmup ? c.logAverage : logE)));
    } else {
      c.retained++;
      if (accepted) c.accepted++;
      if (!Number.isFinite(delta) || Math.abs(delta) > 1000) c.divergences++;
    }
    return {
      q: c.q.slice(),
      path: proposal.path,
      accepted,
      alpha,
      delta,
      warmup: warming,
      iteration: c.iteration
    };
  }
  // Component labels are a presentation transform, never an ordering constraint
  // inside HMC. One permutation carries every parameter of a component together.
  function mixtureComponents(q) {
    const w = sigmoid(q[4]);
    return [0, 1].map(i => ({
      index: i,
      mu: wrap(q[2 * i]),
      etaMu: q[2 * i],
      etaK: q[2 * i + 1],
      k: softplus(q[2 * i + 1]),
      weight: i ? 1 - w : w,
      etaWeight: i ? -q[4] : q[4]
    })).sort((a, b) => a.mu - b.mu || a.index - b.index);
  }

  function treatmentComponents(q) {
    return [0, 1].map(i => ({
      index: i,
      mu: wrap(q[2 * i]),
      etaMu: q[2 * i],
      etaK: q[2 * i + 1],
      k: softplus(q[2 * i + 1]),
      weight: 1
    }));
  }

  function logadd(a, b) {
    const hi = Math.max(a, b);
    return hi + Math.log1p(Math.exp(Math.min(a, b) - hi));
  }
  // A von Mises prior is defined on the circular mean, not a Normal
  // prior on its predictor. Normalize the auxiliary winding representation
  // so its pushforward is exactly von Mises while the real-line lift is proper.
  // For winding scale 2π, higher Fourier terms are < 1e-33.
  function circularMeanPrior(angle, location = 0, concentration = 0) {
    const a = 2 * Math.exp(-.5 * TAU * TAU),
      normalizer = 1 + a * Math.cos(angle),
      derivative = -a * Math.sin(angle);
    return {
      logp: concentration * Math.cos(angle - location) - (concentration ? bessel(concentration).logI0 :
        0) - .5 * (angle / TAU) ** 2 - Math.log(TAU * Math.sqrt(2 * Math.PI)) - Math.log(normalizer),
      grad: -concentration * Math.sin(angle - location) - angle / (TAU * TAU) - derivative / normalizer
    };
  }

  function mixture(data) {
    const n = data.length;

    function target(q) {
      const [a, e, b, f, z] = q, k = softplus(e), l = softplus(f), B = bessel(k), D = bessel(l), w =
        sigmoid(z);
      const lw = -softplus(-z),
        lv = -softplus(z);
      let logp = 0;
      const grad = [0, 0, 0, 0, 0];
      for (const x of data) {
        const u = k * Math.cos(x - a) - B.logI0 + lw,
          v = l * Math.cos(x - b) - D.logI0 + lv,
          h = logadd(u, v),
          r = Math.exp(u - h),
          s = Math.exp(v - h);
        logp += h;
        grad[0] += r * k * Math.sin(x - a);
        grad[1] += r * (Math.cos(x - a) - B.A);
        grad[2] += s * l * Math.sin(x - b);
        grad[3] += s * (Math.cos(x - b) - D.A);
        grad[4] += r - w;
      }
      const pa = circularMeanPrior(a),
        pb = circularMeanPrior(b);
      // The weight prior is Normal(0,1) on the shared logit itself.
      // No sigmoid Jacobian belongs in this unconstrained target.
      logp += pa.logp + pb.logp - .5 * (e / 3) ** 2 - .5 * (f / 3) ** 2 - .5 * z * z;
      grad[0] += pa.grad;
      grad[2] += pb.grad;
      grad[1] = sigmoid(e) * grad[1] - e / 9;
      grad[3] = sigmoid(f) * grad[3] - f / 9;
      grad[4] -= z;
      return {
        logp,
        grad
      };
    }
    // Multi-start circular k-means initializes dispersed chains; it does not
    // assign observations during inference, which marginalizes membership.
    let best = null;
    for (let start = 0; start < 12; start++) {
      let means = [-Math.PI + start * TAU / 12, start * TAU / 12],
        sets;
      for (let j = 0; j < 15; j++) {
        sets = [
          [],
          []
        ];
        data.forEach(x => sets[Math.cos(x - means[0]) >= Math.cos(x - means[1]) ? 0 : 1].push(x));
        if (sets.some(s => s.length < 2)) break;
        means = sets.map(s => model(s).mean);
      }
      if (sets.some(s => s.length < 2)) continue;
      const parts = sets.map(s => model(s)),
        q = [means[0], invsoftplus(Math.min(12, parts[0].khat)), means[1], invsoftplus(Math.min(12,
          parts[1].khat)), Math.log(sets[0].length / sets[1].length)];
      const score = target(q).logp;
      if (!best || score > best.score) best = {
        q,
        parts,
        score
      };
    }
    if (!best) {
      const m = model(data);
      best = {
        q: [m.mean, 1, wrap(m.mean + Math.PI), 1, 0],
        parts: [m, m]
      };
    }

    function initialize(random, index) {
      const q = best.q.slice(),
        mass = [];
      for (let i = 0; i < 2; i++) {
        const m = best.parts[i],
          k = softplus(q[2 * i + 1]),
          A = bessel(k).A;
        mass.push(Math.max(1, m.n * k * m.R), Math.max(.2, m.n * Math.max(1e-4, 1 - A * A - A / k) *
          sigmoid(q[2 * i + 1]) ** 2 + 1 / 9));
        q[2 * i] += [-TAU, 0, TAU, 0][index % 4] + .15 * normal(random);
        q[2 * i + 1] += .25 * normal(random);
      }
      mass.push(Math.max(1, n * .2));
      q[4] += .2 * normal(random);
      if (index % 2) {
        const a = q[0],
          e = q[1];
        q[0] = q[2];
        q[1] = q[3];
        q[2] = a;
        q[3] = e;
        q[4] = -q[4];
      }
      return {
        q,
        mass
      };
    }
    return {
      kind: 'mixture',
      n,
      target,
      initialize,
      components: mixtureComponents,
      khat: Math.max(...best.parts.map(m => m.khat))
    };
  }

  function treatments(data, groups) {
    const parts = [0, 1].map(i => model(data.filter((_, j) => groups[j] === i)));
    return {
      kind: 'treatments',
      n: data.length,
      parts,
      khat: Math.max(...parts.map(m => m.khat)),
      components: treatmentComponents,
      target(q) {
        const a = parts[0].target(q.slice(0, 2)),
          b = parts[1].target(q.slice(2, 4));
        return {
          logp: a.logp + b.logp,
          grad: [...a.grad, ...b.grad]
        };
      },
      initialize(random, index) {
        const cs = parts.map(m => chain(m, Math.floor(random() * 4294967296), index));
        return {
          q: cs.flatMap(c => c.q),
          mass: cs.flatMap(c => c.mass)
        };
      }
    };
  }

  function example(data, kind = 'single', groups = []) {
    if (kind === 'mixture') return mixture(data);
    if (kind === 'uniform') return vonMisesUniform(data);
    if (kind === 'treatments') return treatments(data, groups);
    return {
      ...model(data),
      kind: 'single',
      components: q => [{
        index: 0,
        mu: wrap(q[0]),
        etaMu: q[0],
        etaK: q[1],
        k: softplus(q[1]),
        weight: 1
      }]
    };
  }

  function vonMisesUniform(data) {
    function target(q) {
      const [m, e, z] = q, k = softplus(e), B = bessel(k), w = sigmoid(z), lw = -softplus(-z), lu = -
        softplus(z);
      let logp = 0,
        gm = 0,
        gk = 0,
        gz = 0;
      for (const angle of data) {
        // Both components have the same 1/(2π) factor; marginalize Z.
        const v = k * Math.cos(angle - m) - B.logI0 + lw,
          h = logadd(v, lu),
          r = Math.exp(v - h);
        logp += h;
        gm += r * k * Math.sin(angle - m);
        gk += r * (Math.cos(angle - m) - B.A);
        gz += r - w;
      }
      const prior = circularMeanPrior(m),
        center = Math.log(2);
      // Normal(log(2),1) on the logit favors the von Mises component.
      logp += prior.logp - .5 * (e / 3) ** 2 - .5 * (z - center) ** 2;
      return {
        logp,
        grad: [gm + prior.grad, sigmoid(e) * gk - e / 9, gz - (z - center)]
      };
    }
    const fit = model(data);
    let best = null;
    for (const k of [1, 3, 5, 8, 12])
      for (const w of [.3, .5, .7, .9]) {
        const q = [fit.mean, invsoftplus(k), Math.log(w / (1 - w))],
          score = target(q).logp;
        if (!best || score > best.score) best = {
          q,
          score
        };
      }
    return {
      kind: 'uniform',
      n: data.length,
      khat: softplus(best.q[1]),
      target,
      components(q) {
        const w = sigmoid(q[2]);
        return [{
          index: 0,
          mu: wrap(q[0]),
          etaMu: q[0],
          etaK: q[1],
          k: softplus(q[1]),
          weight: w,
          etaWeight: q[2]
        }, {
          index: 1,
          uniform: true,
          weight: 1 - w
        }];
      },
      initialize(random, index) {
        const q = best.q.slice(),
          k = softplus(q[1]),
          A = bessel(k).A,
          w = sigmoid(q[2]);
        const mass = [Math.max(1, data.length * w * k * A), Math.max(.2, data.length * w * Math.max(
          1e-4, 1 - A * A - A / k) * sigmoid(q[1]) ** 2 + 1 / 9), Math.max(1, data.length * w * (
          1 - w))];
        q[0] += [-TAU, 0, TAU, 0][index % 4] + .25 * normal(random);
        q[1] += .3 * normal(random);
        q[2] += .25 * normal(random);
        return {
          q,
          mass
        };
      }
    };
  }
  const api = {
    TAU,
    wrap,
    softplus,
    sigmoid,
    invsoftplus,
    rng,
    normal,
    bessel,
    vonMises,
    model,
    leapfrog,
    chain,
    step,
    mixtureComponents,
    treatmentComponents,
    mixture,
    treatments,
    vonMisesUniform,
    example,
    circularMeanPrior
  };
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.CircularInference = api;
})(typeof self !== 'undefined' ? self : globalThis);
