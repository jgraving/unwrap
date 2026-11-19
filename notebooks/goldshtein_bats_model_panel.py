"""Prototype a route map and fitted heading check for bat figure panel b."""
from __future__ import annotations

import argparse
from pathlib import Path
from types import SimpleNamespace

import goldshtein_bats_layout_review as composite
import goldshtein_bats_posthoc as posthoc
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import xarray as xr

from goldshtein_bats_layout_data import _goal_window_mass, load_layout_data
from unwrap.link.numpy import circular_modulo, softplus


PAGE_W, PAGE_H = 438, 269
BIN_EDGES = np.arange(0, 181, 15)
OUTPUT_DIR = composite.source.ROOT / 'assets' / 'goldshtein-bats-figure-ideas'
MODEL_STYLES = (('von Mises', '#68737a', '--'), ('von Mises–Uniform', '#312441', '-'))
OBS_COLOR = '#d8dcde'


def _baseline_parameters(posterior, tags, treatment_indices, levels):
    """Evaluate the pinned Bambi baseline predictors for the observed bats."""

    def predictor(prefix):
        intercept = posthoc._stack_samples_first(posterior[prefix + 'Intercept'])
        effects = posterior[prefix + 'treatment']
        effects = effects.sel({effects.dims[-1]: np.asarray(levels[1:]).astype(str)})
        offsets = posthoc._stack_samples_first(effects)
        offsets = np.column_stack([np.zeros_like(intercept), offsets])
        random_effects = posterior[prefix + '1|tag'].sel({'tag__factor_dim': tags.astype(str)})
        return (intercept[:, None] + offsets[:, treatment_indices]
                + posthoc._stack_samples_first(random_effects))

    return circular_modulo(predictor('')), softplus(predictor('kappa_'))


def _fitted_bin_mass(mu, kappa, weight, step_shares):
    """Integrate absolute error bins and average over the modeled steps."""
    edges = np.deg2rad(BIN_EDGES)
    by_draw = np.zeros((len(mu), len(edges) - 1))
    for i, share in enumerate(step_shares):
        cdf = _goal_window_mass(edges, mu[:, i, None], kappa[:, i, None],
                                weight[:, i, None], posthoc.KAPPA_UNIFORM)
        by_draw += share * np.diff(cdf, axis=1)
    np.testing.assert_allclose(by_draw.sum(axis=1), 1, atol=2e-12)
    lower, upper = posthoc.hdi_along_samples(by_draw, prob=.95)
    return by_draw.mean(axis=0), lower, upper


def heading_check(data_root, layout):
    """Compare observed and fitted bin probabilities for the same 4,942 steps."""
    rows = layout['tracks'].loc[layout['evidence'].index]
    tags, first, counts = np.unique(rows['tag'], return_index=True, return_counts=True)
    levels = layout['levels']
    treatment_indices = np.searchsorted(levels, rows['treatment'].to_numpy()[first])
    shares = counts / counts.sum()
    baseline_path = data_root / 'fit_cache' / 'bats_vm_v1.nc'
    mixture_path = data_root / 'fit_cache' / 'bats_vmu_v1.nc'

    with xr.open_dataset(baseline_path, group='observed_data', engine='h5netcdf') as observed:
        angles = observed['HeadingRadians'].values
    with xr.open_dataset(mixture_path, group='observed_data', engine='h5netcdf') as observed:
        np.testing.assert_allclose(angles, observed['Heading_obs'].values, atol=1e-12)
    errors = np.deg2rad(rows['step_bearing_deg'] - rows['bearing_to_dest_deg'])
    np.testing.assert_allclose(np.abs(circular_modulo(errors)), np.abs(angles), atol=1e-12)
    observed_mass = np.histogram(np.rad2deg(np.abs(angles)), BIN_EDGES)[0] / len(angles)

    with xr.open_dataset(baseline_path, group='posterior', engine='h5netcdf') as posterior:
        mu, kappa = _baseline_parameters(posterior, tags, treatment_indices, levels)
        baseline = _fitted_bin_mass(mu, kappa, np.broadcast_to(1, mu.shape), shares)
    with xr.open_dataset(mixture_path, group='posterior', engine='h5netcdf') as posterior:
        fit = SimpleNamespace(posterior=posterior)
        posthoc.restore_per_obs_deterministics(fit, np.arange(len(tags)), treatment_indices)
        mu, kappa, weight = posthoc.extract_simulation_arrays(fit, n_draws=None)
        mixture = _fitted_bin_mass(mu, kappa, weight, shares)
    return observed_mass, (baseline, mixture)


def _draw_routes(ax, layout):
    tracks = layout['tracks']
    colors = dict(zip(layout['levels'], composite.PALETTE, strict=True))
    lake = layout['lake']
    polygons = [lake] if lake.geom_type == 'Polygon' else list(lake.geoms)
    features = np.vstack([tracks[['x_m', 'y_m']], layout['sites'], layout['roost'],
                          np.reshape(lake.bounds, (2, 2))])
    xlim = features[:, 0].min() - 225, features[:, 0].max() + 225
    ylim = features[:, 1].min() - 300, features[:, 1].max() + 300
    ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], facecolor=composite.MAP_FACE)
    ax.set_aspect('equal')
    for polygon in polygons:
        ax.fill(*polygon.exterior.xy, fc='#d9e9ed', ec='#a3bcc3', lw=.5, zorder=1)
    for treatment, group in tracks.groupby('treatment', sort=True, observed=True):
        segments, gaps = [], []
        for _, flight in group.groupby('tag', sort=False, observed=True):
            xy = flight[['x_m', 'y_m']].to_numpy()
            edges = np.stack((xy[:-1], xy[1:]), axis=1)
            gap = flight['gap_before'].to_numpy()[1:] == 1
            segments.extend(edges[~gap])
            gaps.extend(edges[gap])
        color = colors[treatment]
        ax.add_collection(LineCollection(segments, colors=color, lw=.64, alpha=.62,
                                         capstyle='round', zorder=3))
        ax.add_collection(LineCollection(gaps, colors=color, lw=.5, alpha=.4,
                                         linestyles=(0, (2, 2)), zorder=2))
    for point in layout['sites']:
        ax.plot(*point, marker='X', ms=4.5, mfc=composite.INK, mec='white', mew=.65, zorder=4)
    ax.plot(*layout['roost'], marker='*', ms=6.4, mfc=composite.INK, mec='white', mew=.65, zorder=4)
    x, y = xlim[1] - 430, ylim[0] + 450
    ax.plot([x - 1000, x], [y, y], color=composite.INK, lw=1.1)
    ax.annotate('1 km', (x - 500, y), xytext=(0, 3), textcoords='offset points',
                ha='center', fontsize=7.5)
    north = ylim[1] - 1320
    ax.annotate('', (x, north + 625), (x, north),
                arrowprops={'arrowstyle': '-|>', 'lw': .8, 'color': composite.INK,
                            'shrinkA': 0, 'shrinkB': 0})
    ax.annotate('N', (x, north + 625), xytext=(0, 4), textcoords='offset points',
                ha='center', fontsize=7.5)


def draw_model_panel(layout, check):
    """Place the observed flights beside the two fitted heading distributions."""
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72), facecolor='white')

    def axes_at(x, y, w, h):
        return fig.add_axes([x / PAGE_W, 1 - (y + h) / PAGE_H, w / PAGE_W, h / PAGE_H])

    def page_text(x, y, text, **kwargs):
        return fig.text(x / PAGE_W, 1 - y / PAGE_H, text, va='center', **kwargs)

    page_text(12, 20, 'b.', fontsize=16, weight='bold')
    page_text(34, 20, 'Homing trajectories and model fit', fontsize=12, weight='bold')
    page_text(63, 43, 'Routes by treatment', ha='center', fontsize=10, weight='bold')
    _draw_routes(axes_at(12, 52, 102, 184), layout)
    for label, color, x, y in zip(composite.LABELS, composite.PALETTE,
                                  (18, 69, 18, 69), (246, 246, 258, 258), strict=True):
        page_text(x, y, label, color=color, fontsize=8)

    observed, fits = check
    fit = axes_at(167, 70, 254, 111)
    residual = axes_at(167, 196, 254, 43)
    page_text(294, 43, 'Observed and fitted headings', ha='center', fontsize=10, weight='bold')
    handles = [Patch(facecolor=OBS_COLOR, edgecolor='none', label='Observed')]
    handles.extend(Line2D([], [], color=color, ls=style, label=name)
                   for name, color, style in MODEL_STYLES)
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(294 / PAGE_W, 1 - 57 / PAGE_H),
               ncols=3, frameon=False, fontsize=8, columnspacing=1.2, handlelength=1.7)

    fit.bar(BIN_EDGES[:-1], observed, width=15, align='edge', color=OBS_COLOR, linewidth=0)
    residual.axhline(0, color='#9a9a9a', lw=.6)
    for (mean, lo, hi), (_, color, style) in zip(fits, MODEL_STYLES, strict=True):
        for ax, reference in ((fit, np.zeros_like(observed)), (residual, observed)):
            ax.stairs(mean - reference, BIN_EDGES, baseline=None, color=color, ls=style, lw=1.2)
            ax.fill_between(BIN_EDGES, np.r_[lo - reference, (lo - reference)[-1]],
                            np.r_[hi - reference, (hi - reference)[-1]], step='post',
                            color=color, alpha=.12, linewidth=0)
    for ax in (fit, residual):
        ax.set(xlim=(0, 180), xticks=(0, 45, 90, 135, 180))
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(length=2.8, width=.65, labelsize=8)
    fit.set(ylim=(0, None), ylabel='Probability per 15° bin')
    fit.tick_params(axis='x', bottom=False, labelbottom=False)
    residual.set(xlabel='Absolute error from goal bearing (°)', ylabel='Fitted −\nobserved')
    return fig


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path,
                        default=composite.source.ROOT / '.cache/goldshtein-bats-figure' /
                        composite.source.DATA_REVISION)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    parser.add_argument('--font-dir', type=Path)
    parser.add_argument('--font-family', default='Arial')
    args = parser.parse_args(argv)
    root = composite.source.prepare_inputs(args.data_root, download=args.download)
    layout = load_layout_data(root)
    check = heading_check(root, layout)
    composite.source.resolve_fonts(args.font_dir, args.font_family)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {'font.family': args.font_family, 'font.size': 9, 'pdf.fonttype': 42,
             'text.color': composite.INK, 'axes.labelcolor': composite.INK,
             'axes.edgecolor': composite.INK, 'axes.linewidth': .65}
    with mpl.rc_context(style):
        composite.save_figure(draw_model_panel(layout, check), args.output_dir, 'bat_model_panel')


if __name__ == '__main__':
    main()
