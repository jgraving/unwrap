/* Canvas views share one worker-backed sample stream; rendering never samples. */
(() => {
  'use strict';
  const C = window.CircularInference,
    $ = id => document.getElementById(id),
    TAU = C.TAU,
    pi = Math.PI;
  const ids = ['data', 'predictive', 'raw', 'linked', 'trace-mu', 'trace-k', 'trace-weight', 'effects',
    'concentration-effects'
  ];
  const canvases = Object.fromEntries(ids.map(id => [id, $('circ-' + id)]));
  if (!canvases.data) return;
  const script = document.currentScript.src,
    worker = new Worker(new URL('worker.js', script));
  const motion = matchMedia('(prefers-reduced-motion:reduce)');
  const state = {
    predictor: 0,
    individual: 0,
    design: null,
    offsets: [],
    kOffsets: [],
    kind: 'single',
    component: 0,
    groups: [],
    inputs: {},
    seed: 42,
    generation: 0,
    data: [],
    model: null,
    chains: [
      [],
      [],
      [],
      []
    ],
    stats: [],
    positions: [],
    last: [],
    iteration: 0,
    running: false,
    ready: false,
    complete: false,
    traces: 'raw',
    stamp: 0,
    dirty: true,
    gridDirty: true,
    predictive: null,
    proposalsAt: 0
  };
  let palette = {},
    grids = {},
    raf = 0,
    lastPaint = 0,
    liveFrame = 0,
    resumeAfterReset = false;
  const conditional = () => state.kind === 'regression' || state.kind === 'hierarchical';
  const sampleContours = () => conditional() || state.kind === 'mixture' || state.kind === 'uniform';
  const controls = ['n', 'mu', 'kappa', 'mu2', 'kappa2', 'weight'].map(id => $('circ-' + id));
  const fmt = x => Math.abs(x) >= 10 ? x.toFixed(0) : x.toFixed(1),
    degrees = x => Math.round(x * 180 / pi) + '°';

  function colors() {
    const s = getComputedStyle(document.documentElement);
    const get = k => s.getPropertyValue(k).trim();
    palette = {
      ink: get('--ink'),
      muted: get('--muted'),
      line: get('--line'),
      paper: get('--surface'),
      blue: get('--circ-blue'),
      teal: get('--circ-teal'),
      orange: get('--circ-orange'),
      mauve: get('--circ-mauve')
    };
    palette.chains = [palette.blue, palette.teal, palette.orange, palette.mauve];
  }
  colors();
  const chainVisible = () => true;

  function notify(text) {
    $('circ-status').textContent = text;
  }
  const icons = {
    "play": "<svg viewBox=\"0 0 24 24\" width=\"22\" height=\"22\" aria-hidden=\"true\"><path d=\"M9 5.5 19 12 9 18.5Z\" fill=\"currentColor\"/></svg>",
    "pause": "<svg viewBox=\"0 0 24 24\" width=\"22\" height=\"22\" aria-hidden=\"true\"><path d=\"M8 6v12M16 6v12\" fill=\"none\" stroke=\"currentColor\" stroke-width=\"3\" stroke-linecap=\"butt\"/></svg>"
  };

  function playbackIcon() {
    const b = document.querySelector('[data-action=run]'),
      label = state.complete ? 'Sampling complete' : state.running ? 'Pause sampling' :
      'Play sampling';
    b.innerHTML = state.running ? icons.pause : icons.play;
    b.setAttribute('aria-label', label);
    b.title = label;
    b.setAttribute('aria-pressed', String(state.running));
  }

  function stop() {
    state.running = false;
    worker.postMessage({
      type: 'pause'
    });
    playbackIcon();
  }

  function invalidate() {
    state.dirty = true;
    if (!raf) raf = requestAnimationFrame(paint);
  }

  function reset(newData) {
    stop();
    state.generation++;
    state.ready = false;
    state.complete = false;
    state.chains = [
      [],
      [],
      [],
      []
    ];
    state.iteration = 0;
    state.stats = [];
    state.positions = [];
    state.last = [];
    state.predictive = null;
    if (newData) {
      if (state.kind === 'mixture') {
        const ordered = [{
          mu: C.wrap(+controls[1].value * pi / 180) * 180 / pi,
          k: controls[2].value,
          w: +controls[5].value
        }, {
          mu: C.wrap(+controls[3].value * pi / 180) * 180 / pi,
          k: controls[4].value,
          w: 1 - controls[5].value
        }].sort((a, b) => a.mu - b.mu);
        controls[1].value = Math.round(ordered[0].mu);
        controls[2].value = ordered[0].k;
        controls[3].value = Math.round(ordered[1].mu);
        controls[4].value = ordered[1].k;
        controls[5].value = ordered[0].w;
        updateValues();
      }
      const random = C.rng(state.seed),
        n = +controls[0].value,
        mu = +controls[1].value * pi / 180,
        k = +controls[2].value,
        mu2 = +controls[3].value * pi / 180,
        k2 = +controls[4].value,
        weight = +controls[5].value;
      if (conditional()) {
        const generated = C.generatePredictorData(state.kind, {
          n,
          mu,
          slope: mu2,
          k,
          k2,
          spread: +$('circ-individual-sd').value * pi / 180,
          kSpread: +$('circ-individual-k-sd').value,
          varyConcentration: $('circ-vary-concentration').checked,
          concentrationSlope: +$('circ-concentration-slope').value,
          seed: state.seed
        });
        state.data = generated.data;
        state.groups = generated.groups;
        state.design = generated.design;
        state.offsets = generated.offsets;
        state.kOffsets = generated.kOffsets;
      } else {
        state.design = null;
        state.groups = [];
        state.data = Array.from({
          length: n
        }, (_, i) => {
          const group = state.kind === 'treatments' ? (i < Math.ceil(n / 2) ? 0 : 1) : (state
            .kind === 'mixture' || state.kind === 'uniform') ? (random() < weight ? 0 : 1) : 0;
          state.groups.push(group);
          if (state.kind === 'uniform' && group) return C.wrap(TAU * (random() - .5));
          return C.vonMises(group ? mu2 : mu, group ? k2 : k, random);
        });
      }

    }
    state.model = C.example(state.data, state.kind, state.groups, state.design);
    const k = state.model.khat;
    state.kmax = Math.max(5, Math.min(35, k + Math.max(4, 3 * k * Math.sqrt(2 / (state.data.length / (
      state.kind === 'single' ? 1 : 2))))));
    grids = {};
    state.gridDirty = true;
    $('circ-seed').textContent = 'Dataset ' + state.seed;
    document.querySelector('[data-action=run]').disabled = true;
    notify('Preparing chains…');
    worker.postMessage({
      type: 'init',
      generation: state.generation,
      seed: state.seed + 100003,
      data: state.data,
      kind: state.kind,
      groups: state.groups,
      design: state.design
    });
    playbackIcon();
    updateStats();
    invalidate();
  }

  function updateStats() {
    const n = state.stats.reduce((s, x) => s + x.retained, 0),
      accepted = state.stats.reduce((s, x) => s + x.accepted, 0),
      divs = state.stats.reduce((s, x) => s + x.divergences, 0);
    $('circ-draws').textContent = n.toLocaleString();
    $('circ-acceptance').textContent = n ? Math.round(100 * accepted / n) + '%' : '—';
    $('circ-divergences').textContent = divs;
    $('circ-phase').textContent = state.complete ? 'Complete' : state.iteration === 0 ? 'Ready' : state
      .iteration <= 120 ? 'Warmup' : 'Sampling';
    $('circ-mixing').textContent = divs ?
      'Divergences suggest that the sampler may have missed part of the posterior. Interpret these estimates cautiously and compare the chains on the constrained scale.' :
      'Different unconstrained values can represent the same direction on a circle. Compare constrained traces when assessing whether the chains agree.';
  }
  worker.onmessage = ({
    data: d
  }) => {
    if (d.generation !== state.generation) return;
    if (d.type === 'ready') {
      worker.postMessage({
        type: 'speed',
        value: +$('circ-speed').value
      });
      state.ready = true;
      document.querySelector('[data-action=run]').disabled = false;
      notify('Ready');
      state.positions = d.positions;
      if (resumeAfterReset) {
        resumeAfterReset = false;
        state.running = true;
        playbackIcon();
        worker.postMessage({
          type: 'run'
        });
        notify('Warmup');
      }
      invalidate();
      return;
    }
    if (d.type !== 'draws') return;
    for (const row of d.rows) {
      row.forEach((s, i) => {
        state.chains[i].push(s);
        state.positions[i] = s.q;
      });
      state.last = row;
      state.iteration = row[0].iteration;
    }
    state.stats = d.stats;
    state.complete = d.complete;
    state.gridDirty = true;
    state.predictive = null;
    state.proposalsAt = performance.now();
    if (d.complete) {
      stop();
      document.querySelector('[data-action=run]').disabled = true;
    }
    if (state.complete) notify('Complete · 1,000 draws / chain');
    else if (state.iteration === 120 || state.iteration === 121 || state.iteration % 20 < 12) notify(
      state.iteration <= 120 ? 'Warmup · ' + state.iteration + ' / 120' : 'Sampling · ' + Math.max(
        0, state.iteration - 120) + ' / 1,000');
    updateStats();
    invalidate();
  };
  worker.onerror = () => {
    stop();
    state.ready = false;
    document.querySelector('[data-action=run]').disabled = true;
    notify('The sampler could not load. Reload this example to try again.');
  };
  document.querySelectorAll('[data-action]').forEach(b => b.addEventListener('click', () => {
    const action = b.dataset.action;
    if (action === 'generate') {
      state.seed++;
      reset(true);
    } else if (action === 'restart') reset(false);
    else if (action === 'run' && state.ready && !state.complete) {
      if (state.running) {
        stop();
        notify('Paused');
      } else {
        state.running = true;
        playbackIcon();
        worker.postMessage({
          type: 'run'
        });
        notify(state.iteration < 120 ? 'Warmup' : 'Sampling');
      }
    }
  }));

  function updateValues() {
    ['n', 'mu', 'kappa', 'mu2', 'kappa2', 'weight'].forEach((name, i) => {
      $('circ-' + name + '-value').textContent = name === 'n' ? controls[i].value : name
        .startsWith('mu') ? controls[i].value + '°' : name === 'weight' ? Math.round(+controls[i]
          .value * 100) + '%' : (+controls[i].value).toFixed(1);
    });
  }

  function updateTraceTitles() {
    const linked = state.traces === 'linked';
    $('circ-mean-trace-title').innerHTML = linked ? 'Mean angle μ' : 'Mean predictor η<sub>μ</sub>';
    $('circ-k-trace-title').innerHTML = linked ? 'Concentration κ' :
      'Concentration predictor η<sub>κ</sub>';
    $('circ-weight-trace-title').innerHTML = linked ? 'Mixture weight w' :
      'Weight predictor η<sub>w</sub>';
  }

  function updateExample() {
    const multi = state.kind === 'mixture' || state.kind === 'treatments',
      mixture = state.kind === 'mixture',
      uniform = state.kind === 'uniform',
      weighted = mixture || uniform;
    document.documentElement.dataset.circularExample = state.kind;
    document.querySelectorAll('[data-multi]').forEach(e => {
      e.hidden = !multi;
    });
    document.querySelectorAll('[data-mixture]').forEach(e => {
      e.hidden = !weighted;
    });
    $('circ-component-control').hidden = !multi;
    $('circ-group-key').hidden = state.kind !== 'treatments';
    $('circ-weight-label').textContent = uniform ? 'von Mises weight' : 'Mixing weight';
    $('circ-component-1-key').textContent = uniform ? 'w × von Mises' : 'w × component 1';
    $('circ-component-2-key').textContent = uniform ? '(1 − w) × uniform' : '(1 − w) × component 2';
    $('circ-component-label').textContent = mixture ? 'Component' : 'Treatment';
    const options = $('circ-component').options;
    for (let i = 0; i < 2; i++) options[i].textContent = mixture ? 'Component ' + (i + 1) :
      'Treatment ' + (i ? 'B' : 'A');
    $('circ-mu-label').textContent = multi ? (mixture ? 'Mean 1' : 'Mean A') : 'Mean angle';
    $('circ-kappa-label').textContent = multi ? (mixture ? 'Concentration 1' : 'Concentration A') :
      'Concentration';
    $('circ-mu2-label').textContent = mixture ? 'Mean 2' : 'Mean B';
    $('circ-kappa2-label').textContent = mixture ? 'Concentration 2' : 'Concentration B';
    $('circ-example-note').textContent = uniform ?
      'Infer a directional signal mixed with uniform background.' : mixture ?
      'Infer two directions, two concentrations and one mixing weight.' : state.kind === 'treatments' ?
      'Compare two groups, each with its own mean and concentration.' :
      'Infer a mean direction and concentration.';
    $('circ-model-note').textContent = uniform ?
      'This example separates a directional signal from a uniform background. The sampler estimates the signal’s direction, concentration and share of the observations, while allowing uncertainty about all three.' :
      mixture ?
      'This example allows two directional components, each with its own mean direction and concentration. One mixing weight describes their relative contributions. Components are kept in a consistent angular order so their parameters and curves stay aligned.' :
      state.kind === 'treatments' ?
      'This example compares two known treatment groups. Each group has its own mean direction and concentration. Both groups are fitted together; the treatment selector changes which group’s posterior geometry and chains you inspect.' :
      'This example estimates one population’s mean direction and concentration. Concentration describes how tightly observations cluster around the mean direction.';
    $('circ-prior-note').textContent = weighted ?
      'The model combines prior assumptions about plausible directions, concentrations and component weights with the evidence in the observations. Circular priors are used for mixture directions; the signal–background example modestly favors a directional signal before seeing the data.' :
      'The model combines prior assumptions about plausible parameter values with the evidence in the observations. Prior values in an unconstrained space are mapped into valid circular-model parameters before they are compared with the data.';
    $('circ-contour-label').textContent = weighted ? 'Sample contours' : 'Posterior contours';
    $('circ-raw').setAttribute('aria-label', uniform ?
      'Marginal posterior of the von Mises mean and concentration predictors, with sample density, estimated contours and MCMC paths.' :
      mixture ?
      'Marginal posterior for the selected sorted component: unconstrained mean and concentration predictors, sampled density, estimated sample contours and MCMC paths.' :
      'Joint posterior in unconstrained mean and concentration predictors, showing exact posterior contours, sampled density and MCMC paths.'
      );
    $('circ-data').setAttribute('aria-label', uniform ?
      'Unlabeled angular observations and the generating von Mises–uniform mixture density.' :
      mixture ? 'Unlabeled observations and the generating two-component von Mises mixture.' : state
      .kind === 'treatments' ?
      'Observed directions for treatments A and B with their generating densities.' :
      'Generated angular observations around a circle with their generating von Mises density.');
    $('circ-predictive').setAttribute('aria-label', uniform ?
      'Posterior von Mises–uniform mixture densities and a pointwise 90 percent credible band.' :
      mixture ? 'Posterior mixture densities and a pointwise 90 percent credible band.' : state
      .kind === 'treatments' ?
      'Separate posterior densities and pointwise 90 percent bands for treatments A and B.' :
      'Von Mises posterior densities and a pointwise 90 percent credible band.');
    if (conditional()) {
      const hierarchical = state.kind === 'hierarchical';
      document.querySelectorAll('[data-multi]').forEach(e => {
        e.hidden = true;
      });
      $('circ-mu2').closest('.circ-control').hidden = false;
      $('circ-component-control').hidden = !hierarchical;
      $('circ-group-key').hidden = !hierarchical;
      $('circ-mu-label').textContent = hierarchical ? 'Mean A' : 'Mean at x = 0';
      $('circ-mu2-label').textContent = hierarchical ? 'Treatment effect' : 'Change per unit x';
      $('circ-kappa-label').textContent = hierarchical ? 'Group concentration A' :
        'Concentration at x = 0';
      $('circ-kappa2').closest('.circ-control').hidden = !hierarchical;
      $('circ-kappa2-label').textContent = 'Group concentration B';
      $('circ-example-note').textContent = hierarchical ?
        'Both mean direction and concentration vary by treatment and include partially pooled individual differences.' :
        'Mean direction changes with a continuous predictor; concentration can stay constant or vary with x.';
      $('circ-model-note').textContent = hierarchical ?
        'Each of four individuals contributes observations in both treatments. A shared treatment effect shifts mean direction, while individual offsets come from one learned distribution. Concentration also depends on treatment, with individual offsets drawn from a second learned distribution. Both sets of individual effects are partially pooled across the four individuals.' :
        'A linear predictor sets mean direction: μ(x) = wrap(α + βx). Concentration is constant by default. With varying concentration, κ(x) = softplus(γ₀ + γ₁x). Both predictors are estimated jointly from all observations.';
      $('circ-prior-note').textContent = hierarchical ?
        'The mean intercept has a Normal(0, 2π) prior and the treatment effect a Normal(0, 1.5) prior, in radians. Standardized individual offsets have independent Normal(0, 1) priors; their shared spread has a half-Normal(0, 0.8) prior. Concentration follows κᵢₜ = softplus(γₜ + σₖvᵢ): treatment predictors γₜ have independent Normal(0, 3) priors, standardized individual offsets vᵢ have Normal(0, 1) priors, and their shared spread σₖ has a half-Normal(0, 1) prior. Mean-direction and concentration offsets are independent a priori. The sampler uses a non-centered parameterization, including the spread transformation’s Jacobian.' :
        'The intercept has a Normal(0, 2π) prior and the slope a Normal(0, 1.5) prior, in radians per predictor unit. The concentration intercept has a Normal(0, 3) prior; its optional slope has a Normal(0, 1.5) prior. Softplus maps their linear predictor to positive concentration. Priors are defined on these unconstrained predictors; displayed wrapped-angle densities are posterior projections.';
      $('circ-contour-label').textContent = 'Sample contours';
      $('circ-effects-title').textContent = hierarchical ? 'Treatments and individual differences' :
        'Mean direction across the predictor';
      $('circ-concentration-effects-title').textContent = hierarchical ?
        'Hierarchical concentration by treatment' : 'Concentration across the predictor';
      $('circ-concentration-note').textContent =
        'Colored estimates show partially pooled individuals, with the selected individual emphasized. The bold line shows the group baseline; bars show 90% credible intervals. The dashed line generated the group baseline.';
      $('circ-effects-note').textContent = hierarchical ?
        'Thin colored lines show individual means; the bold line shows treatment means before individual offsets. Bars show 90% credible intervals.' :
        'Dots show observations; dashed curves generated the data. Solid curves and shaded bands show posterior means and pointwise 90% credible intervals. The direction band wraps at ±180°.';
      $('circ-raw').setAttribute('aria-label',
        'Sample-based marginal posterior for the selected mean direction and concentration, in unconstrained coordinates.'
        );
      $('circ-data').setAttribute('aria-label', hierarchical ?
        'Observed directions for the selected individual in both treatments.' :
        'Observed directions near the selected predictor value.');
      $('circ-predictive').setAttribute('aria-label', hierarchical ?
        'Conditional posterior densities for the selected individual in both treatments.' :
        'Conditional posterior density at the selected predictor value.');
    }
    document.querySelectorAll('[data-conditional]').forEach(e => {
      e.hidden = !conditional();
    });
    document.querySelectorAll('[data-regression]').forEach(e => {
      e.hidden = state.kind !== 'regression';
    });
    document.querySelectorAll('[data-hierarchical]').forEach(e => {
      e.hidden = state.kind !== 'hierarchical';
    });
    $('circ-concentration-slope-control').hidden = state.kind !== 'regression' || !$(
      'circ-vary-concentration').checked;
    updateConditionalSelection();
    updateTraceTitles();
    updateValues();
  }

  function liveData() {
    updateValues();
    if (liveFrame) return;
    resumeAfterReset = resumeAfterReset || state.running || state.complete;
    liveFrame = requestAnimationFrame(() => {
      liveFrame = 0;
      reset(true);
    });
  }
  controls.forEach(input => {
    input.addEventListener('input', liveData);
  });
  $('circ-example').addEventListener('change', e => {
    if (liveFrame) {
      cancelAnimationFrame(liveFrame);
      liveFrame = 0;
    }
    resumeAfterReset = false;
    state.inputs[state.kind] = controls.slice(1).map(c => c.value);
    state.kind = e.target.value;
    state.component = 0;
    $('circ-component').value = '0';
    if (conditional() && !state.inputs[state.kind]) controls[0].value = '80';
    const values = state.inputs[state.kind] || (state.kind === 'regression' ? [145, 4, 65, 4, .5] :
      state.kind === 'hierarchical' ? [-25, 3, 65, 8, .5] : state.kind === 'uniform' ? [60, 5, 0,
        1, .5
      ] : state.kind === 'mixture' ? [-65, 3, 75, 8, .65] : state.kind === 'treatments' ? [-40,
        2, 60, 6, .5
      ] : [155, 1, -145, 1, .5]);
    controls.slice(1).forEach((c, i) => {
      c.value = values[i];
    });
    updateExample();
    reset(true);
  });
  $('circ-component').addEventListener('change', e => {
    state.component = +e.target.value;
    grids = {};
    state.gridDirty = true;
    invalidate();
  });
  $('circ-speed').addEventListener('change', e => worker.postMessage({
    type: 'speed',
    value: +e.target.value
  }));
  $('circ-traces').addEventListener('change', e => {
    state.traces = e.target.value;
    updateTraceTitles();
    invalidate();
  });
  document.addEventListener('visibilitychange', () => {
    if (document.hidden && state.running) {
      stop();
      notify('Paused');
    }
  });
  addEventListener('pagehide', () => {
    worker.terminate();
    cancelAnimationFrame(raf);
  });
  const observer = new ResizeObserver(invalidate);
  Object.values(canvases).forEach(c => observer.observe(c));
  new MutationObserver(() => {
    colors();
    invalidate();
  }).observe(document.documentElement, {
    attributes: true,
    attributeFilter: ['data-theme', 'data-motion']
  });
  motion.addEventListener('change', invalidate);

  function setup(id) {
    const canvas = canvases[id],
      r = canvas.getBoundingClientRect(),
      dpr = Math.min(devicePixelRatio || 1, 2.5),
      w = Math.max(1, r.width),
      h = Math.max(1, r.height);
    if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
      canvas.width = Math.round(w * dpr);
      canvas.height = Math.round(h * dpr);
    }
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.font = '11px ' + getComputedStyle(document.documentElement).getPropertyValue('--font').trim();
    ctx.lineCap = 'round';
    ctx.lineJoin = 'round';
    return {
      ctx,
      w,
      h
    };
  }

  function text(ctx, s, x, y, color = palette.muted, align = 'center') {
    ctx.save();
    ctx.fillStyle = color;
    if (!/η[μκw]/.test(s)) {
      ctx.textAlign = align;
      ctx.fillText(s, x, y);
      ctx.restore();
      return;
    }
    const font = ctx.font,
      size = parseFloat(font) || 11,
      parts = [],
      pattern = /η([μκw])/g;
    let cursor = 0,
      m;
    while ((m = pattern.exec(s))) {
      if (m.index > cursor) parts.push({
        text: s.slice(cursor, m.index),
        sub: false
      });
      parts.push({
        text: 'η',
        sub: false
      }, {
        text: m[1],
        sub: true
      });
      cursor = m.index + m[0].length;
    }
    if (cursor < s.length) parts.push({
      text: s.slice(cursor),
      sub: false
    });
    for (const part of parts) {
      ctx.font = part.sub ? font.replace(/[\d.]+px/, Math.round(size * .75) + 'px') : font;
      part.width = ctx.measureText(part.text).width;
    }
    const width = parts.reduce((sum, p) => sum + p.width, 0);
    let left = x - (align === 'center' ? width / 2 : align === 'right' ? width : 0);
    ctx.textAlign = 'left';
    for (const part of parts) {
      ctx.font = part.sub ? font.replace(/[\d.]+px/, Math.round(size * .75) + 'px') : font;
      ctx.fillText(part.text, left, y + (part.sub ? size * .3 : 0));
      left += part.width;
    }
    ctx.restore();
  }

  function path(ctx, points, color, width = 1, alpha = 1, fill = false) {
    if (!points.length) return;
    ctx.save();
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    points.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]));
    if (fill) {
      ctx.closePath();
      ctx.fillStyle = color;
      ctx.fill();
    } else {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.stroke();
    }
    ctx.restore();
  }

  function dot(ctx, x, y, color, r = 3, alpha = 1) {
    ctx.save();
    ctx.globalAlpha = alpha;
    ctx.beginPath();
    ctx.arc(x, y, r, 0, TAU);
    ctx.fillStyle = color;
    ctx.fill();
    ctx.restore();
  }

  function axes(view, xlim, ylim, xlabel, ylabel, xticks, yticks) {
    const {
      ctx,
      w,
      h
    } = view, left = 50, right = 16, top = 12, bottom = 40;
    const X = x => left + (x - xlim[0]) / (xlim[1] - xlim[0]) * (w - left - right),
      Y = y => h - bottom - (y - ylim[0]) / (ylim[1] - ylim[0]) * (h - bottom - top);
    ctx.strokeStyle = palette.line;
    ctx.lineWidth = 1;
    for (const y of yticks) {
      path(ctx, [
        [left, Y(y)],
        [w - right, Y(y)]
      ], palette.line, .7, .55);
      text(ctx, fmt(y), left - 7, Y(y) + 4, palette.muted, 'right');
    }
    for (const [x, label] of xticks) {
      path(ctx, [
        [X(x), h - bottom],
        [X(x), h - bottom + 4]
      ], palette.muted, .7);
      text(ctx, label, X(x), h - bottom + 17);
    }
    path(ctx, [
      [left, top],
      [left, h - bottom],
      [w - right, h - bottom]
    ], palette.muted, .8);
    text(ctx, xlabel, (left + w - right) / 2, h - 3);
    ctx.save();
    ctx.translate(12, (top + h - bottom) / 2);
    ctx.rotate(-pi / 2);
    text(ctx, ylabel, 0, 0);
    ctx.restore();
    return {
      X,
      Y,
      left,
      right,
      top,
      bottom,
      w,
      h,
      ctx
    };
  }

  function components(q) {
    return state.model.components(q, {
      x: state.predictor,
      individual: state.individual
    });
  }

  function selected(q) {
    return components(q)[state.component];
  }

  function project(q) {
    const c = selected(q);
    return [c.etaMu, c.etaK];
  }

  function retained(full = false) {
    return state.chains.flatMap(xs => xs.filter(x => !x.warmup).map(x => full ? x.q : project(x.q)));
  }
  const pdf = (angle, mu, k) => Math.exp(k * Math.cos(angle - mu) - C.bessel(k).logI0) / TAU;

  function drawCircle(id, predictive) {
    const {
      ctx,
      w,
      h
    } = setup(id), cx = w / 2, cy = h / 2 + 1, R = Math.min(w, h) * .34;
    const polar = (a, r) => [cx + Math.sin(a) * r, cy - Math.cos(a) * r];
    const angles = Array.from({
        length: 181
      }, (_, i) => -pi + TAU * i / 180),
      treatment = state.kind === 'treatments' || state.kind === 'hierarchical',
      mixture = state.kind === 'mixture' || state.kind === 'uniform';
    const truthParts = [{
      mu: +controls[1].value * pi / 180,
      k: +controls[2].value,
      weight: mixture ? +controls[5].value : 1
    }, {
      mu: +controls[3].value * pi / 180,
      k: +controls[4].value,
      weight: mixture ? 1 - controls[5].value : 1
    }];
    if (conditional()) {
      const offset = state.kind === 'hierarchical' ? state.offsets[state.individual] : 0,
        base = +controls[1].value * pi / 180,
        beta = +controls[3].value * pi / 180;
      const kOffset = state.kind === 'hierarchical' ? state.kOffsets[state.individual] : 0;
      truthParts[0] = {
        mu: base + offset + (state.kind === 'regression' ? beta * state.predictor : 0),
        k: state.kind === 'regression' ? generatingK(state.predictor) : C.softplus(C.invsoftplus(+
          controls[2].value) + kOffset),
        weight: 1
      };
      truthParts[1] = {
        mu: base + beta + offset,
        k: state.kind === 'hierarchical' ? C.softplus(C.invsoftplus(+controls[4].value) + kOffset) :
          +controls[2].value,
        weight: 1
      };
    }
    const displayData = state.data.map((angle, index) => ({
      angle,
      index
    })).filter(({
        index
      }) => state.kind === 'regression' ? Math.abs(state.design.x[index] - state.predictor) <= .25 :
      state.kind === 'hierarchical' ? state.design.individual[index] === state.individual : true);
    const truth = angles.map(a => pdf(a, truthParts[0].mu, truthParts[0].k));
    const truth2 = angles.map(a => state.kind === 'uniform' ? 1 / TAU : pdf(a, truthParts[1].mu,
      truthParts[1].k));
    const weightedTruth = mixture ? [truth.map(v => v * truthParts[0].weight), truth2.map(v => v *
      truthParts[1].weight)] : null;
    const truths = treatment ? [truth, truth2] : [mixture ? weightedTruth[0].map((v, i) => v +
      weightedTruth[1][i]) : truth];
    // Physical polar radius: r = 1 + density, with one shared fixed radial limit.
    const densityLimit = Math.ceil(Math.max(...truths.flat(), conditional() ? pdf(0, 0, Math.max(+
        controls[2].value, state.kind === 'hierarchical' ? +controls[4].value : 0)) : .05) * 1.12 *
      20) / 20;
    const radialLimit = 1 + densityLimit,
      unit = R / radialLimit,
      base = unit;
    const radius = v => unit * (1 + v);
    for (const a of [0, pi / 2, pi, 3 * pi / 2]) path(ctx, [polar(a, base), polar(a, R + 12)], palette
      .line, .7);
    ctx.strokeStyle = palette.muted;
    ctx.lineWidth = .9;
    ctx.beginPath();
    ctx.arc(cx, cy, base, 0, TAU);
    ctx.stroke();
    text(ctx, '0°', cx, cy - R - 17);
    text(ctx, '90°', cx + R + 16, cy + 4);
    text(ctx, '±180°', cx, cy + R + 24);
    text(ctx, '−90°', cx - R - 16, cy + 4);
    // Density ticks start at zero on the unit circle; never label physical r as density.
    const axis = pi / 4;
    path(ctx, [polar(axis, base), polar(axis, radius(densityLimit))], palette.line, .65, .7);
    for (const value of [0, densityLimit / 2, densityLimit]) {
      const [x, y] = polar(axis, radius(value));
      path(ctx, [
        [x - 2, y - 2],
        [x + 2, y + 2]
      ], palette.muted, .65, .7);
      text(ctx, String(Number(value.toFixed(2))), x + 5, y - 4, palette.muted, 'left');
    }
    if (predictive && !state.predictive) {
      const samples = retained(true);
      // Sample the full span of all four chains even during the first few post-warmup updates.
      const count = Math.min(samples.length, 100);
      const batch = Array.from({length: count}, (_, index) => {
        const position = Math.floor((index + 0.5) * samples.length / count);
        return components(samples[position]);
      });
      if (batch.length) {
        state.predictive = Array.from({
          length: treatment ? 2 : 1
        }, (_, group) => {
          // Weight every draw before averaging; E[w f] is not E[w] E[f].
          // The component permutation has already been applied to each full draw.
          const contributions = mixture ? batch.map(cs => cs.map(c => {
            const logI0 = c.uniform ? 0 : C.bessel(c.k).logI0;
            return angles.map(a => c.weight * (c.uniform ? 1 / TAU : Math.exp(c.k * Math.cos(
              a - c.mu) - logI0) / TAU));
          })) : null;
          const curves = contributions ? contributions.map(cs => angles.map((_, i) => cs[0][i] +
            cs[1][i])) : batch.map(cs => angles.map(a => treatment ? pdf(a, cs[group].mu, cs[
            group].k) : pdf(a, cs[0].mu, cs[0].k)));
          const componentMeans = contributions ? [0, 1].map(j => angles.map((_, i) => contributions
            .reduce((sum, cs) => sum + cs[j][i], 0) / contributions.length)) : null;
          const mean = [],
            lo = [],
            hi = [];
          angles.forEach((a, i) => {
            const xs = curves.map(c => c[i]).sort((a, b) => a - b);
            mean.push(xs.reduce((s, x) => s + x, 0) / xs.length);
            lo.push(xs[Math.floor(.05 * (xs.length - 1))]);
            hi.push(xs[Math.floor(.95 * (xs.length - 1))]);
          });
          return {
            curves,
            mean,
            lo,
            hi,
            componentMeans
          };
        });
      }
    }
    if (predictive && !state.predictive) return;
    // Fit the generating PDF with modest headroom, never the evolving posterior band.
    // Both circles therefore share one linear density scale throughout sampling.
    const P = predictive ? state.predictive : [];
    const coords = values => values.map((v, i) => polar(angles[i], radius(v)));
    // A fixed polar rlim clips outlying posterior curves instead of moving the truth.
    ctx.save();
    ctx.beginPath();
    ctx.arc(cx, cy, R, 0, TAU);
    ctx.clip();
    if (predictive) {
      P.forEach((p, group) => {
        const color = mixture ? palette.ink : group ? palette.teal : palette.blue;
        path(ctx, [...coords(p.hi), ...coords(p.lo).reverse()], color, 0, mixture ? .08 : .13,
          true);
        p.curves.filter((_, i) => i % 8 === 0).forEach(c => path(ctx, coords(c), color, .8,
          mixture ? .055 : .1));
        if (mixture) p.componentMeans.forEach((c, j) => {
          const componentColor = j ? palette.teal : palette.blue;
          path(ctx, [...angles.map(a => polar(a, base)), ...coords(c).reverse()],
            componentColor, 0, .07, true);
          path(ctx, coords(c), componentColor, 1.6, .9);
        });
        path(ctx, coords(p.mean), color, mixture ? 2.5 : 2.2);
      });
    } else {
      if (mixture) {
        weightedTruth.forEach((c, j) => {
          const color = j ? palette.teal : palette.blue;
          path(ctx, [...angles.map(a => polar(a, base)), ...coords(c).reverse()], color, 0, .07,
            true);
          ctx.save();
          ctx.setLineDash([4, 3]);
          path(ctx, coords(c), color, 1.6, .9);
          ctx.restore();
        });
      } else truths.forEach((t, group) => path(ctx, [...angles.map(a => polar(a, base)), ...coords(t)
        .reverse()
      ], group ? palette.teal : palette.blue, 0, .06, true));
      const totals = new Map(),
        counts = new Map(),
        binFor = (a, group) => ((Math.round(a * 18 / pi) + 36) % 36) + 36 * group;
      displayData.forEach(({
        angle: a,
        index: i
      }) => {
        const bin = binFor(a, treatment ? state.groups[i] : 0);
        totals.set(bin, (totals.get(bin) || 0) + 1);
      });
      displayData.forEach(({
        angle: a,
        index: i
      }) => {
        const group = treatment ? state.groups[i] : 0,
          bin = binFor(a, group),
          j = counts.get(bin) || 0,
          n = totals.get(bin);
        counts.set(bin, j + 1);
        // Two treatments use disjoint bounded annuli, preserving exact angles.
        const outer = Math.max(1, base - 5) * (treatment && group ? .79 : 1),
          inner = outer * (treatment ? .85 : .72),
          r = n > 1 ? outer - (outer - inner) * j / (n - 1) : outer;
        dot(ctx, ...polar(a, r), mixture ? palette.muted : group ? palette.teal : palette.blue,
          Math.min(2.5, base * .055), .8);
      });
      if (state.kind === 'single') {
        const end = polar(state.model.mean, base * .9);
        path(ctx, [
          [cx, cy], end
        ], palette.blue, 2);
        const v = [end[0] - cx, end[1] - cy],
          length = Math.hypot(...v),
          u = [v[0] / length, v[1] / length];
        path(ctx, [end, [end[0] - 8 * u[0] + 4 * u[1], end[1] - 8 * u[1] - 4 * u[0]],
          [end[0] - 8 * u[0] - 4 * u[1], end[1] - 8 * u[1] + 4 * u[0]]
        ], palette.blue, 0, 1, true);
      }
    }
    ctx.save();
    ctx.setLineDash([4, 4]);
    truths.forEach((t, group) => path(ctx, coords(t), mixture && !predictive ? palette.ink :
      treatment ? (group ? palette.teal : palette.blue) : palette.muted, mixture && !predictive ?
      2 : 1.2, mixture && !predictive ? .9 : .7));
    ctx.restore();
    ctx.restore();
  }
  // Derive view limits from retained draws, never from warmup excursions.
  function sampleLimits(values, fallback, lower = -Infinity, upper = Infinity) {
    const finite = values.filter(Number.isFinite);
    if (!finite.length) return fallback.slice();
    let lo = Math.min(...finite),
      hi = Math.max(...finite),
      span = Math.max(hi - lo, .1);
    const padding = span * .12;
    lo = Math.max(lower, lo - padding);
    hi = Math.min(upper, hi + padding);
    const rough = (hi - lo) / 4,
      power = 10 ** Math.floor(Math.log10(rough)),
      step = ([1, 2, 5, 10].find(n => n * power >= rough) || 10) * power;
    return [Math.max(lower, Math.floor(lo / step) * step), Math.min(upper, Math.ceil(hi / step) *
      step)];
  }
  const rangeTicks = range => Array.from({
    length: 5
  }, (_, i) => range[0] + (range[1] - range[0]) * i / 4);

  function sampleDomain(linked) {
    const draws = retained(),
      x = draws.map(q => q[0]),
      y = draws.map(q => linked ? C.softplus(q[1]) : q[1]);
    return {
      xlim: linked ? [-pi, pi] : sampleLimits(x, [-3 * pi, 3 * pi]),
      ylim: sampleLimits(y, linked ? [0, state.kmax] : [-5, C.invsoftplus(state.kmax)], linked ? 0 : -
        Infinity)
    };
  }

  function posteriorGrid(linked, domain = sampleDomain(linked)) {
    const nx = 100,
      ny = 68,
      {
        xlim,
        ylim
      } = domain,
      values = [];
    let max = -Infinity;
    for (let y = 0; y < ny; y++)
      for (let x = 0; x < nx; x++) {
        const a = xlim[0] + (x + .5) / nx * (xlim[1] - xlim[0]),
          b = ylim[0] + (y + .5) / ny * (ylim[1] - ylim[0]);
        const model = state.kind === 'treatments' ? state.model.parts[state.component] : state.model;
        const v = sampleContours() ? 0 : linked ? model.constrained(a, b) : model.target([a, b]).logp;
        values.push(v);
        max = Math.max(max, v);
      }
    const weights = values.map(x => Math.exp(x - max)),
      sorted = [...weights].sort((a, b) => b - a),
      total = sorted.reduce((a, b) => a + b, 0),
      levels = [];
    let sum = 0,
      index = 0;
    for (const coverage of [.5, .8, .95]) {
      while (sum < coverage * total && index < sorted.length) sum += sorted[index++];
      levels.push(sorted[Math.max(0, index - 1)]);
    }
    return {
      nx,
      ny,
      xlim,
      ylim,
      weights,
      levels
    };
  }
  // Marching squares on numerically evaluated density; no smoothing of axes.
  function contours(ctx, g, A) {
    const {
      nx,
      ny,
      xlim,
      ylim,
      weights,
      levels
    } = g;
    const point = (x, y) => [A.X(xlim[0] + (x + .5) / nx * (xlim[1] - xlim[0])), A.Y(ylim[0] + (y +
      .5) / ny * (ylim[1] - ylim[0]))];
    levels.forEach((level, index) => {
      ctx.beginPath();
      for (let y = 0; y < ny - 1; y++)
        for (let x = 0; x < nx - 1; x++) {
          const vs = [weights[y * nx + x], weights[y * nx + x + 1], weights[(y + 1) * nx + x + 1],
              weights[(y + 1) * nx + x]
            ],
            ps = [
              [x, y],
              [x + 1, y],
              [x + 1, y + 1],
              [x, y + 1]
            ],
            cuts = [];
          for (let e = 0; e < 4; e++) {
            const f = (e + 1) % 4;
            if ((vs[e] >= level) === (vs[f] >= level)) continue;
            const t = (level - vs[e]) / (vs[f] - vs[e]);
            cuts.push(point(ps[e][0] + t * (ps[f][0] - ps[e][0]), ps[e][1] + t * (ps[f][1] - ps[e][
              1
            ])));
          }
          for (let j = 0; j + 1 < cuts.length; j += 2) {
            ctx.moveTo(...cuts[j]);
            ctx.lineTo(...cuts[j + 1]);
          }
        }
      ctx.strokeStyle = palette.muted;
      ctx.globalAlpha = .5 - index * .1;
      ctx.lineWidth = index === 0 ? 1.2 : .8;
      ctx.stroke();
      ctx.globalAlpha = 1;
    });
  }

  function histogram(g, linked) {
    const counts = new Float32Array(g.nx * g.ny),
      {
        xlim,
        ylim,
        nx,
        ny
      } = g;
    retained().forEach(q => {
      const a = linked ? C.wrap(q[0]) : q[0],
        b = linked ? C.softplus(q[1]) : q[1];
      const x = Math.floor((a - xlim[0]) / (xlim[1] - xlim[0]) * nx),
        y = Math.floor((b - ylim[0]) / (ylim[1] - ylim[0]) * ny);
      if (x < 0 || x >= nx || y < 0 || y >= ny) return;
      for (let dy = -2; dy <= 2; dy++)
        for (let dx = -2; dx <= 2; dx++) {
          const yy = y + dy;
          let xx = x + dx;
          if (linked) xx = (xx + nx) % nx;
          if (xx < 0 || xx >= nx || yy < 0 || yy >= ny) continue;
          counts[yy * nx + xx] += Math.exp(-(dx * dx + dy * dy) / 2);
        }
    });
    g.counts = counts;
    g.maxCount = Math.max(...counts, .01);
    if (sampleContours()) {
      g.weights = Array.from(counts);
      const sorted = Array.from(counts).sort((a, b) => b - a),
        total = sorted.reduce((s, x) => s + x, 0);
      let sum = 0,
        index = 0;
      g.levels = total ? [.5, .8, .95].map(p => {
        while (sum < p * total && index < sorted.length) sum += sorted[index++];
        return sorted[Math.max(0, index - 1)];
      }) : [];
    }
  }

  function drawJoint(id, linked, now) {
    const view = setup(id),
      domain = sampleDomain(linked),
      cached = grids[id];
    if (!cached || [...domain.xlim, ...domain.ylim].some((v, i) => v !== [...cached.xlim, ...cached
        .ylim
      ][i])) grids[id] = posteriorGrid(linked, domain);
    const g = grids[id];
    if (state.gridDirty || !g.counts) histogram(g, linked);
    const {
      ctx
    } = view, {
      xlim,
      ylim,
      nx,
      ny
    } = g;
    const A = axes(view, xlim, ylim, linked ? 'Mean angle μ' : 'Mean predictor ημ', linked ?
      'Concentration κ' : 'Concentration predictor ηκ', linked ? [
        [-pi, '−180°'],
        [-pi / 2, '−90°'],
        [0, '0°'],
        [pi / 2, '90°'],
        [pi, '180°']
      ] : rangeTicks(xlim).map(x => [x, degrees(x)]), rangeTicks(ylim));
    ctx.save();
    ctx.beginPath();
    ctx.rect(A.left, A.top, A.w - A.left - A.right, A.h - A.top - A.bottom);
    ctx.clip();
    const cw = (A.w - A.left - A.right) / nx,
      ch = (A.h - A.top - A.bottom) / ny;
    for (let y = 0; y < ny; y++)
      for (let x = 0; x < nx; x++) {
        const d = g.counts[y * nx + x] / g.maxCount;
        if (d < .005) continue;
        ctx.globalAlpha = Math.min(.78, Math.sqrt(d) * .8);
        ctx.fillStyle = d > .45 ? palette.teal : palette.blue;
        ctx.fillRect(A.left + x * cw, A.h - A.bottom - (y + 1) * ch, cw + .5, ch + .5);
      }
    ctx.globalAlpha = 1;
    contours(ctx, g, A);
    const transform = q => {
      const c = selected(q);
      return [linked ? c.mu : c.etaMu, linked ? c.k : c.etaK, c.index];
    };

    function trajectory(qs, color, width, alpha, dashed = false) {
      ctx.save();
      if (dashed) ctx.setLineDash([3, 3]);
      let segment = [],
        previous = null;
      for (const q of qs) {
        const t = transform(q);
        if (previous && ((linked && Math.abs(t[0] - previous[0]) > pi) || (state.kind === 'mixture' &&
            t[2] !== previous[2]))) {
          path(ctx, segment, color, width, alpha);
          segment = [];
        }
        segment.push([A.X(t[0]), A.Y(t[1])]);
        previous = t;
      }
      path(ctx, segment, color, width, alpha);
      ctx.restore();
    }
    state.chains.forEach((xs, i) => {
      if (!chainVisible(i)) return;
      trajectory(xs.slice(-30).map(s => s.q), palette.chains[i], .9, .4);
      const last = state.last[i];
      if (last) {
        const fraction = motion.matches ? 1 : Math.min(1, (now - state.proposalsAt) / 60);
        const length = Math.max(2, Math.round(last.path.length * fraction));
        trajectory(last.path.slice(0, length), palette.chains[i], 1.5, last.accepted ? .9 : .35, !
          last.accepted);
      }
      const q = state.positions[i];
      if (q) {
        const t = transform(q);
        dot(ctx, A.X(t[0]), A.Y(t[1]), palette.paper, 4.7);
        dot(ctx, A.X(t[0]), A.Y(t[1]), palette.chains[i], 3);
      }
    });
    ctx.restore();
  }

  function drawTrace(id, index) {
    const view = setup(id),
      linked = state.traces === 'linked',
      xmax = Math.max(180, Math.ceil(state.iteration / 100) * 100);
    const value = q => {
      const c = selected(q);
      return index === 2 ? (linked ? c.weight : c.etaWeight) : index === 1 ? (linked ? c.k : c
        .etaK) : (linked ? c.mu : c.etaMu);
    };
    const values = state.chains.flatMap(xs => xs.filter(s => !s.warmup).map(s => value(s.q)));
    let ylim;
    if (index === 0) ylim = linked ? [-pi, pi] : [-3 * pi, 3 * pi];
    else if (index === 1) ylim = linked ? [0, state.kmax] : [-5, C.invsoftplus(state.kmax)];
    else ylim = linked ? [0, 1] : [-4, 4];
    if (!(linked && index === 0)) ylim = sampleLimits(values, ylim, linked && index ? 0 : -Infinity,
      linked && index === 2 ? 1 : Infinity);
    const label = index === 2 ? (linked ? 'w' : 'ηw') : index === 1 ? (linked ? 'κ' : 'ηκ') : (linked ?
      'μ (degrees)' : 'ημ (degrees)');
    const A = axes(view, [0, xmax], ylim, 'Iteration', label, [
      [0, '0'],
      [120, '120'],
      [xmax, '' + xmax]
    ], index ? [ylim[0], (ylim[0] + ylim[1]) / 2, ylim[1]] : []);
    const {
      ctx
    } = view;
    if (index === 0)
      for (const y of linked ? [-pi, 0, pi] : rangeTicks(ylim)) text(ctx, '' + Math.round(y * 180 /
        pi), A.left - 6, A.Y(y) + 4, palette.muted, 'right');
    ctx.save();
    ctx.beginPath();
    ctx.rect(A.left, A.top, A.w - A.left - A.right, A.h - A.top - A.bottom);
    ctx.clip();
    ctx.fillStyle = palette.muted;
    ctx.globalAlpha = .07;
    ctx.fillRect(A.X(0), A.top, A.X(120) - A.X(0), view.h - A.top - A.bottom);
    ctx.globalAlpha = 1;
    state.chains.forEach((xs, i) => {
      let points = [],
        last = null,
        lastComponent = null;
      for (const s of xs) {
        const v = value(s.q),
          component = selected(s.q).index;
        if ((linked && index === 0 && last !== null && Math.abs(v - last) > pi) || (state.kind ===
            'mixture' && lastComponent !== null && component !== lastComponent)) {
          path(ctx, points, palette.chains[i], .85, .75);
          points = [];
        }
        points.push([A.X(s.iteration), A.Y(v)]);
        last = v;
        lastComponent = component;
      }
      path(ctx, points, palette.chains[i], .85, .75);
    });
    ctx.restore();
  }

  function generatingK(x) {
    return $('circ-vary-concentration').checked ? C.softplus(C.invsoftplus(+controls[2].value) + +$(
      'circ-concentration-slope').value * x) : +controls[2].value;
  }
  $('circ-vary-concentration').addEventListener('change', () => {
    updateExample();
    liveData();
  });
  $('circ-concentration-slope').addEventListener('input', event => {
    $('circ-concentration-slope-value').textContent = (+event.target.value).toFixed(1);
    liveData();
  });

  function updateConditionalSelection() {
    if (!$('circ-conditional-selection')) return;
    $('circ-predictor-value').textContent = state.predictor.toFixed(2);
    const where = 'x = ' + state.predictor.toFixed(2);
    $('circ-data-heading').textContent = 'Observed directions';
    $('circ-predictive-heading').textContent = state.kind === 'regression' ? 'Posterior' :
      'Posterior on the circle';
    for (const id of ['circ-data-condition', 'circ-predictive-condition']) {
      $(id).hidden = state.kind !== 'regression';
      $(id).textContent = where;
    }
    $('circ-conditional-selection').textContent = state.kind === 'regression' ?
      'Polar plots and posterior geometry show ' + (state.predictor === 0 ? 'the intercept, x = 0' :
        'x = ' + state.predictor.toFixed(2)) +
      '. Observed headings come from nearby predictor values.' : state.kind === 'hierarchical' ?
      'Circles and posterior geometry follow Individual ' + (state.individual + 1) +
      ' in both treatments.' : '';
  }

  function reproject() {
    state.predictive = null;
    grids = {};
    state.gridDirty = true;
    updateConditionalSelection();
    invalidate();
  }
  $('circ-predictor').addEventListener('input', event => {
    state.predictor = +event.target.value;
    reproject();
  });
  $('circ-individual').addEventListener('change', event => {
    state.individual = +event.target.value;
    reproject();
  });
  $('circ-individual-sd').addEventListener('input', event => {
    $('circ-individual-sd-value').textContent = event.target.value + '°';
  });
  $('circ-individual-sd').addEventListener('input', liveData);
  $('circ-individual-k-sd').addEventListener('input', event => {
    $('circ-individual-k-sd-value').textContent = (+event.target.value).toFixed(1);
    liveData();
  });

  function circularSummary(values) {
    if (!values.length) return null;
    const center = Math.atan2(values.reduce((s, v) => s + Math.sin(v), 0), values.reduce((s, v) => s +
      Math.cos(v), 0));
    const residuals = values.map(v => C.wrap(v - center)).sort((a, b) => a - b),
      at = p => residuals[Math.floor(p * (residuals.length - 1))];
    return {
      mean: center,
      low: center + at(.05),
      high: center + at(.95)
    };
  }

  function drawEffects() {
    const view = setup('effects'),
      hierarchical = state.kind === 'hierarchical',
      A = axes(view, hierarchical ? [-.35, 1.35] : [-1, 1], [-pi, pi], hierarchical ? 'Treatment' :
        'Predictor x', 'Mean direction', hierarchical ? [
          [0, 'A'],
          [1, 'B']
        ] : [
          [-1, '−1'],
          [-.5, '−0.5'],
          [0, '0'],
          [.5, '0.5'],
          [1, '1']
        ], []),
      ctx = view.ctx;
    for (const y of [-pi, -pi / 2, 0, pi / 2, pi]) text(ctx, degrees(y), A.left - 7, A.Y(y) + 4,
      palette.muted, 'right');
    const draws = retained(true),
      stride = Math.max(1, Math.ceil(draws.length / 160)),
      batch = draws.filter((_, i) => i % stride === 0).slice(-160);

    function line(values, color, width = 2, alpha = 1) {
      let segment = [];
      for (const v of values) {
        if (segment.length && Math.abs(v[1] - segment[segment.length - 1][1]) > pi) {
          path(ctx, segment.map(p => [A.X(p[0]), A.Y(p[1])]), color, width, alpha);
          segment = [];
        }
        segment.push(v);
      }
      path(ctx, segment.map(p => [A.X(p[0]), A.Y(p[1])]), color, width, alpha);
    }

    function interval(x, summary, color, width = 1.2) {
      if (!summary) return;
      const lo = summary.low,
        hi = summary.high;
      for (let j = -1; j <= 1; j++) {
        const a = Math.max(-pi, lo + j * TAU),
          b = Math.min(pi, hi + j * TAU);
        if (b >= a) path(ctx, [
          [A.X(x), A.Y(a)],
          [A.X(x), A.Y(b)]
        ], color, width, .5);
      }
    }
    ctx.save();
    ctx.beginPath();
    ctx.rect(A.left, A.top, view.w - A.left - A.right, view.h - A.top - A.bottom);
    ctx.clip();
    if (!hierarchical) {
      state.data.forEach((angle, i) => dot(ctx, A.X(state.design.x[i]), A.Y(angle), palette.muted, 2,
        .5));
      const grid = Array.from({
          length: 41
        }, (_, i) => -1 + i / 20),
        truth = grid.map(x => [x, C.wrap(+controls[1].value * pi / 180 + +controls[3].value * pi /
          180 * x)]);
      ctx.save();
      ctx.setLineDash([5, 4]);
      line(truth, palette.muted, 1.4, .85);
      ctx.restore();
      const fitted = grid.map(x => ({
        x,
        s: circularSummary(batch.map(q => q[0] + q[1] * x))
      })).filter(v => v.s);
      for (let i = 1; i < fitted.length; i++) {
        const a = fitted[i - 1],
          b = fitted[i],
          shift = TAU * Math.round((a.s.mean - b.s.mean) / TAU);
        for (let turn = -2; turn <= 2; turn++) path(ctx, [
          [A.X(a.x), A.Y(a.s.low + turn * TAU)],
          [A.X(b.x), A.Y(b.s.low + shift + turn * TAU)],
          [A.X(b.x), A.Y(b.s.high + shift + turn * TAU)],
          [A.X(a.x), A.Y(a.s.high + turn * TAU)]
        ], palette.blue, 1, .18, true);
      }
      line(fitted.map(v => [v.x, v.s.mean]), palette.blue, 2.5);
      path(ctx, [
        [A.X(state.predictor), A.top],
        [A.X(state.predictor), view.h - A.bottom]
      ], palette.teal, 1, .55);
    } else {
      const personColors = [palette.blue, palette.teal, palette.orange, palette.mauve];
      for (let person = 0; person < 4; person++) {
        const shift = (person - 1.5) * .09,
          points = [];
        for (let t = 0; t < 2; t++) {
          const indices = state.data.map((_, i) => i).filter(i => state.design.individual[i] ===
            person && state.groups[i] === t);
          indices.forEach((i, j) => dot(ctx, A.X(t + shift + (j % 3 - 1) * .012), A.Y(state.data[i]),
            personColors[person], 1.8, .25));
          const summary = circularSummary(batch.map(q => state.model.components(q, {
            individual: person
          })[t].mu));
          if (summary) {
            interval(t + shift, summary, personColors[person]);
            dot(ctx, A.X(t + shift), A.Y(summary.mean), personColors[person], person === state
              .individual ? 4 : 2.8);
            points.push([t + shift, summary.mean]);
          }
        }
        line(points, personColors[person], person === state.individual ? 2 : 1, .75);
      }
      const means = [];
      for (let t = 0; t < 2; t++) {
        const summary = circularSummary(batch.map(q => q[0] + q[1] * t));
        if (summary) {
          interval(t, summary, palette.ink, 3);
          dot(ctx, A.X(t), A.Y(summary.mean), palette.ink, 4.3);
          means.push([t, summary.mean]);
        }
      }
      line(means, palette.ink, 3);
    }
    ctx.restore();
    if (!draws.length) {
      $('circ-effect-summary').textContent = 'Posterior summaries appear after warmup.';
      return;
    }
    const effects = draws.map(q => q[1] * 180 / pi).sort((a, b) => a - b),
      at = p => effects[Math.floor(p * (effects.length - 1))],
      mean = effects.reduce((a, b) => a + b, 0) / effects.length;
    $('circ-effect-summary').textContent = (hierarchical ? 'Treatment effect: ' : 'Slope: ') + mean
      .toFixed(1) + '°' + (hierarchical ? '' : ' per unit x') + ' · 90% interval ' + at(.05).toFixed(
      1) + '° to ' + at(.95).toFixed(1) + '°' + (hierarchical ? ' · Individual spread: ' + (draws.map(
          q => C.softplus(q[3]) * 180 / pi).sort((a, b) => a - b)[Math.floor(draws.length / 2)])
        .toFixed(1) + '° (posterior median)' : '');
  }

  function drawHierarchicalConcentration() {
    const view = setup('concentration-effects'),
      ctx = view.ctx,
      draws = retained(true),
      stride = Math.max(1, Math.ceil(draws.length / 160)),
      batch = draws.filter((_, i) => i % stride === 0).slice(-160);
    const summary = xs => {
      if (!xs.length) return null;
      xs.sort((a, b) => a - b);
      return {
        low: xs[Math.floor(.05 * (xs.length - 1))],
        high: xs[Math.floor(.95 * (xs.length - 1))],
        median: xs[Math.floor(.5 * (xs.length - 1))]
      };
    };
    const summaries = Array.from({
      length: 5
    }, (_, p) => [0, 1].map(t => summary(batch.map(q => state.model.components(q, {
      individual: p - 1
    })[t].k))));
    const ymax = Math.max(2, +controls[2].value, +controls[4].value, ...summaries.flat().map(s => s
      ?.high || 0)) * 1.12;
    const A = axes(view, [-.35, 1.35], [0, ymax], 'Treatment', 'Concentration κ', [
        [0, 'A'],
        [1, 'B']
      ], [0, ymax / 2, ymax]),
      colors = [palette.ink, ...palette.chains];
    ctx.save();
    ctx.setLineDash([5, 4]);
    path(ctx, [
      [A.X(0), A.Y(+controls[2].value)],
      [A.X(1), A.Y(+controls[4].value)]
    ], palette.muted, 1.4, .8);
    ctx.restore();
    for (let p = 1; p <= 4; p++) {
      const shift = (p - 2.5) * .09,
        points = [];
      for (let t = 0; t < 2; t++) {
        const s = summaries[p][t];
        if (!s) continue;
        const x = t + shift;
        path(ctx, [
          [A.X(x), A.Y(s.low)],
          [A.X(x), A.Y(s.high)]
        ], colors[p], 1.5, .6);
        dot(ctx, A.X(x), A.Y(s.median), colors[p], p - 1 === state.individual ? 4 : 2.8);
        points.push([A.X(x), A.Y(s.median)]);
      }
      path(ctx, points, colors[p], p - 1 === state.individual ? 2 : 1, .75);
    }
    const baseline = [];
    for (let t = 0; t < 2; t++) {
      const s = summaries[0][t];
      if (!s) continue;
      path(ctx, [
        [A.X(t), A.Y(s.low)],
        [A.X(t), A.Y(s.high)]
      ], palette.ink, 3, .65);
      dot(ctx, A.X(t), A.Y(s.median), palette.ink, 4);
      baseline.push([A.X(t), A.Y(s.median)]);
    }
    path(ctx, baseline, palette.ink, 3);
    const out = $('circ-concentration-summary');
    if (!draws.length) {
      out.textContent = 'Hierarchical concentration estimates appear after warmup.';
      return;
    }
    const format = s => s.median.toFixed(1) + ' [' + s.low.toFixed(1) + ', ' + s.high.toFixed(1) + ']';
    const group = [0, 1].map(t => summary(draws.map(q => state.model.components(q, {
        individual: -1
      })[t].k))),
      spread = summary(draws.map(q => C.softplus(q[9]))),
      effect = summary(draws.map(q => q[8] - q[2]));
    out.textContent = 'Group concentration A: ' + format(group[0]) + ' · B: ' + format(group[1]) +
      ' · Treatment effect on ηκ: ' + format(effect) + ' · Individual spread on ηκ: ' + format(
      spread) + ' · Medians [90% credible intervals]';
  }

  function drawConcentration() {
    if (state.kind === 'hierarchical') {
      drawHierarchicalConcentration();
      return;
    }
    const view = setup('concentration-effects'),
      ctx = view.ctx,
      draws = retained(true),
      stride = Math.max(1, Math.ceil(draws.length / 160)),
      batch = draws.filter((_, i) => i % stride === 0).slice(-160),
      grid = Array.from({
        length: 41
      }, (_, i) => -1 + i / 20);
    const quantile = (xs, p) => xs[Math.floor(p * (xs.length - 1))];
    const fitted = grid.map(x => {
      const values = batch.map(q => state.model.components(q, {
        x
      })[0].k).sort((a, b) => a - b);
      return {
        x,
        low: values.length ? quantile(values, .05) : null,
        high: values.length ? quantile(values, .95) : null,
        mean: values.reduce((a, b) => a + b, 0) / Math.max(1, values.length)
      };
    });
    const ymax = Math.max(2, ...grid.map(generatingK), ...fitted.map(v => v.high || 0)) * 1.15,
      A = axes(view, [-1, 1], [0, ymax], 'Predictor x', 'Concentration κ', [
        [-1, '−1'],
        [0, '0'],
        [1, '1']
      ], [0, ymax / 2, ymax]);
    const valid = fitted.filter(v => v.low !== null);
    if (valid.length) {
      path(ctx, [...valid.map(v => [A.X(v.x), A.Y(v.low)]), ...valid.slice().reverse().map(v => [A.X(v
        .x), A.Y(v.high)])], palette.teal, 1, .18, true);
      path(ctx, valid.map(v => [A.X(v.x), A.Y(v.mean)]), palette.teal, 2.5);
    }
    ctx.save();
    ctx.setLineDash([5, 4]);
    path(ctx, grid.map(x => [A.X(x), A.Y(generatingK(x))]), palette.muted, 1.4, .85);
    ctx.restore();
    path(ctx, [
      [A.X(state.predictor), A.top],
      [A.X(state.predictor), view.h - A.bottom]
    ], palette.teal, 1, .55);
  }

  function paint(now) {
    raf = 0;
    if (!state.model) return;
    if (state.dirty && (now - lastPaint > 40 || !state.running)) {
      if (conditional()) drawEffects();
      if (conditional()) drawConcentration();
      drawCircle('data', false);
      drawCircle('predictive', true);
      drawJoint('raw', false, now);
      drawJoint('linked', true, now);
      drawTrace('trace-mu', 0);
      drawTrace('trace-k', 1);
      if (state.kind === 'mixture' || state.kind === 'uniform') drawTrace('trace-weight', 2);
      state.gridDirty = false;
      state.dirty = false;
      lastPaint = now;
    }
    if (state.running && now - state.proposalsAt < 65 && !motion.matches) {
      state.dirty = true;
      raf = requestAnimationFrame(paint);
    }
  }
  updateExample();
  reset(true);
  document.fonts?.ready.then(invalidate);
})();
