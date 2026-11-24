importScripts('core.js', 'regression.js');
const C = self.CircularInference;
let timer = null,
  chains = [],
  model = null,
  seed = 1,
  rate = 1,
  generation = 0;

function pause() {
  clearInterval(timer);
  timer = null;
}

function tick(count = rate) {
  const rows = [];
  for (let i = 0; i < count; i++) {
    if (chains[0].iteration >= 1120) {
      pause();
      break;
    }
    rows.push(chains.map(c => C.step(c, model, 120)));
  }
  postMessage({
    type: 'draws',
    generation,
    rows,
    complete: chains[0].iteration >= 1120,
    stats: chains.map(c => ({
      accepted: c.accepted,
      retained: c.retained,
      divergences: c.divergences,
      epsilon: c.epsilon
    }))
  });
}
self.onmessage = ({
  data: d
}) => {
  if (d.type === 'init') {
    pause();
    generation = d.generation;
    seed = d.seed;
    model = C.example(d.data, d.kind, d.groups, d.design);
    chains = Array.from({
      length: 4
    }, (_, i) => C.chain(model, seed + 7919 * (i + 1), i));
    postMessage({
      type: 'ready',
      generation,
      positions: chains.map(c => c.q),
      mass: chains.map(c => c.mass)
    });
  } else if (d.type === 'run') {
    if (!timer) timer = setInterval(tick, 65);
  } else if (d.type === 'pause') pause();
  else if (d.type === 'step') {
    pause();
    tick(1);
  } else if (d.type === 'speed') rate = d.value;
};
