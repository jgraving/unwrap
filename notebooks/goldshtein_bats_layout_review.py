"""Render the bat composite from the pinned tracks and model fits."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import numpy as np
from PIL import Image
import pymupdf

import goldshtein_bats_publication_figure as source


INK = "#222222"
MAP_FACE = "#f5f6f5"
PALETTE = ("#0072B2", "#BC7300", "#008664", "#C95118")
LABELS = ("Control", "Eye cover", "+ Magnet", "+ Anosmia")
PAGE_W, PAGE_H = 646, 489
OUTPUT_DIR = source.ROOT / 'assets' / 'goldshtein-bats-layout-review'


def draw_composite(data, *, photo_path=None):
    """Compose the photo, matched maps, and three fitted summaries."""
    levels = data["levels"]
    colors = dict(zip(levels, PALETTE, strict=True))
    norm = mpl.colors.Normalize(0, 1)
    cmap = mpl.colormaps["viridis"]
    fig = plt.figure(figsize=(PAGE_W / 72, PAGE_H / 72), facecolor="white")

    def axes_at(x, y, w, h, **kwargs):
        return fig.add_axes([x/PAGE_W, 1-(y+h)/PAGE_H, w/PAGE_W, h/PAGE_H], **kwargs)


    def page_text(x, y, text, **kwargs):
        return fig.text(x/PAGE_W, 1-y/PAGE_H, text, va='center', **kwargs)


    def letter(x, y, value, **kwargs):
        page_text(x, y, value + '.', fontsize=16, weight='bold', **kwargs)


    def setup_axes(ax):
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(length=3, width=.65, pad=3)


    # a. Use the original placeholder and retain all four icon rows as vector art.
    photo = axes_at(22, 25, 151, 145)
    with Image.open(photo_path or source.ASSET_DIR/'image_279.jpeg') as image:
        photo.imshow(image.convert('RGB'))
    photo.axis('off')
    if photo_path is None:
        photo.text(.5, .06, 'Photo placeholder', transform=photo.transAxes,
                   ha='center', va='center', color='white', fontsize=9,
                   bbox={'facecolor': '#272727', 'edgecolor': 'none',
                         'alpha': .58, 'pad': 3.0})
    letter(31, 38, 'a', color='white')
    page_text(22, 188, 'Treatment groups', weight='bold', fontsize=10.8)

    asset = source.load_artwork()
    icons = source.panel_b_drawings(asset)
    ix0, iy0, ix1, iy1 = 202.9, 114.76, 304.79, 234.88
    icon_scale = .76
    icon_top = 200
    icon_ax = axes_at(23, icon_top, (ix1-ix0)*icon_scale, (iy1-iy0)*icon_scale)
    icon_ax.set(xlim=(ix0,ix1), ylim=(source.PAGE_H-iy1,source.PAGE_H-iy0))
    icon_ax.axis('off')
    for index, drawing in sorted(icons.items()):
        patch = source.draw_drawing(icon_ax, drawing, zorder=index+10)
        patch.set_linewidth(patch.get_linewidth()*icon_scale)
    for y, label in zip((126.58,158.73,190.90,223.05), LABELS, strict=True):
        page_text(108, icon_top+(y-iy0)*icon_scale, label, fontsize=9.3)

    # b. Four map panels share metric bounds and a single directedness scale.
    letter(190, 24, 'b')
    page_text(410, 24, 'Homing trajectories', weight='bold', fontsize=11.2, ha='center')
    tracks = data['tracks'].copy()
    tracks['directedness'] = tracks.index.map(data['evidence'])
    track_groups = tracks.groupby('tag', sort=False, observed=True)
    assert len(track_groups) == 67
    assert tracks['directedness'].notna().sum() == len(data['evidence']) == 4942
    assert tracks['directedness'].dropna().between(0, 1).all()
    xy_all = tracks[['x_m','y_m']].to_numpy()
    lake = data['lake']
    polygons = [lake] if lake.geom_type == 'Polygon' else list(lake.geoms)
    features = np.vstack([xy_all, data['sites'], data['roost'], np.reshape(lake.bounds, (2,2))])
    xlim = features[:,0].min()-225, features[:,0].max()+225
    ylim = features[:,1].min()-300, features[:,1].max()+300
    map_h = 216
    map_w = map_h*np.ptp(xlim)/np.ptp(ylim)
    centers = (244, 352, 460, 568)
    counts = {}
    for level, label, center in zip(levels, LABELS, centers, strict=True):
        ax = axes_at(center-map_w/2, 54, map_w, map_h)
        ax.set(xlim=xlim, ylim=ylim, xticks=[], yticks=[], facecolor=MAP_FACE)
        ax.set_aspect('equal')
        for spine in ax.spines.values():
            spine.set_color(colors[level])
            spine.set_linewidth(.8)
        page_text(center, 45, label, ha='center', weight='bold', fontsize=9.8)
        for poly in polygons:
            ax.fill(*poly.exterior.xy, fc='#d9e9ed', ec='#a3bcc3', lw=.55, zorder=2)

        unscored, gaps, segments, means = [], [], [], []
        treatment = tracks.loc[tracks['treatment']==level]
        groups = treatment.groupby('tag', sort=False, observed=True)
        counts[level] = len(groups)
        for _, flight in groups:
            xy = flight[['x_m','y_m']].to_numpy()
            indexed_segments = np.stack((xy[:-1], xy[1:]), axis=1)
            after_gap = flight['gap_before'].to_numpy()[1:] == 1
            # The heading at a row is the displacement ending at that row.
            score = flight['directedness'].to_numpy()[1:]
            valid = np.isfinite(score) & ~after_gap
            gaps.extend(indexed_segments[after_gap])
            unscored.extend(indexed_segments[~np.isfinite(score) & ~after_gap])
            segments.extend(indexed_segments[valid])
            means.extend(score[valid])

        ax.add_collection(LineCollection(unscored, colors='#a7aeae', lw=.53, alpha=.64, zorder=3))
        ax.add_collection(LineCollection(gaps, colors='#959b9c', lw=.5, alpha=.65,
                                         linestyles=(0,(2.3,2)), zorder=3))
        colored = LineCollection(segments, cmap=cmap, norm=norm, lw=.76, alpha=.88,
                                 capstyle='round', zorder=4)
        colored.set_array(np.asarray(means))
        ax.add_collection(colored)
        for site in data['sites']:
            ax.plot(*site, marker='X', ms=5.4, mec='white', mew=.8, mfc=INK, zorder=5)
        ax.plot(*data['roost'], marker='*', ms=7.6, mec='white', mew=.8, mfc=INK, zorder=5)
        if level == levels[-1]:
            x = xlim[1]-430
            y = ylim[0]+470
            ax.plot([x-1000,x],[y,y], color=INK, lw=1.45, zorder=6)
            ax.annotate('1 km', (x-500,y), textcoords='offset points', xytext=(0,4),
                        ha='center', fontsize=8)
            arrow_y = ylim[1]-1470
            ax.annotate('', (x,arrow_y+700), (x,arrow_y),
                        arrowprops=dict(arrowstyle='-|>', lw=1.05, color=INK,
                                        shrinkA=0, shrinkB=0), zorder=6)
            ax.annotate('N',(x,arrow_y+700), textcoords='offset points', xytext=(0,5),
                        ha='center', fontsize=8)
    assert counts == data['counts'] == dict(zip(levels, (29,14,15,9), strict=True))

    cax = axes_at(318, 279, 173, 5.5)
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation='horizontal', ticks=[0,.5,1], format='%g')
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.outline.set_linewidth(.5)
    cb.ax.tick_params(labelsize=8, length=2.25, pad=2)
    page_text(502, 281.6, 'Posterior directedness', fontsize=9.5)

    # c. Step-level PSIS-LOO ELPD and one standard error per model.
    plot_bottom = 456
    letter(22, 314, 'c')
    page_text(103, 314, 'Model comparison', ha='center', weight='bold', fontsize=10.7)
    comparison = axes_at(51, 337, 114, plot_bottom-337)
    for y, key in enumerate(('vm','vmu')):
        score = data['scores'][key]
        comparison.errorbar(score['elpd'], y, xerr=score['se'], color=INK, fmt='o',
                            ms=4.1, elinewidth=1, capsize=2.2, capthick=1, zorder=3)
    lo = min(s['elpd']-s['se'] for s in data['scores'].values())-16
    hi = max(s['elpd']+s['se'] for s in data['scores'].values())+16
    comparison.set(xlim=(lo,hi), ylim=(1.63,-.63), yticks=(0,1),
                   yticklabels=('vM','vMU'), xlabel='ELPD')
    comparison.xaxis.set_major_locator(MultipleLocator(100))
    setup_axes(comparison)
    comparison.spines['left'].set_visible(False)
    comparison.tick_params(axis='y', length=0, pad=4)

    # d. Progression plot. Low-support tails are subdued and defined in the caption.
    letter(182, 314, 'd')
    page_text(309, 314, 'Homing progression', ha='center', weight='bold', fontsize=11)
    c = axes_at(224, 337, 178, plot_bottom-337)
    for level, color in colors.items():
        curve = data['curves'][level]
        x = np.abs(np.asarray(curve['x_mid']))/1000
        mean = np.asarray(curve['mean'])
        support = np.asarray(curve['n_tags_by_x']) >= 5
        c.fill_between(x, curve['lo'], curve['hi'], where=support, color=color,
                       alpha=.14, linewidth=0)
        c.plot(x, np.where(support, mean, np.nan), color=color, lw=1.13, zorder=3)
        c.plot(x, np.where(~support, mean, np.nan), color=color, lw=.85, alpha=.27, zorder=2)
    c.axhline(.5, color='#888888', lw=.65, ls=(0,(2,2)), zorder=0)
    c.set(xlim=(6.6,0), ylim=(0,1), xticks=(6,4,2,0), yticks=(0,.5,1),
          xlabel='Path distance remaining (km)', ylabel='Posterior directedness',
          yticklabels=('0','0.5','1'))
    setup_axes(c)

    # e. Full-mixture polar densities and their 0° contrasts against control.
    letter(420, 314, 'e')
    page_text(443, 314, 'Predictive heading density', weight='bold', fontsize=11)
    bearings = data['density_bearings']
    profiles = data['density_profiles']
    zero_radius = 1.0
    radius_max = zero_radius + 1.08 * max(profile['hi'].max() for profile in profiles.values())
    for center, level, label in zip((444, 497, 550, 603), levels, LABELS, strict=True):
        polar = axes_at(center-23.5, 341, 47, 47, projection='polar')
        polar.set_theta_zero_location('N')
        polar.set_theta_direction(-1)
        polar.set(ylim=(0, radius_max), xticks=[], yticks=[])
        polar.spines['polar'].set_visible(False)
        polar.plot(bearings, np.full_like(bearings, zero_radius), color='#aaaaaa', lw=.55)
        polar.plot([0, 0], [zero_radius, radius_max], color='#9a9a9a', lw=.6)
        profile = profiles[level]
        color = colors[level]
        polar.fill_between(bearings, zero_radius + profile['lo'], zero_radius + profile['hi'],
                           color=color, alpha=.26, linewidth=0)
        polar.plot(bearings, zero_radius + profile['mean'], color=color, lw=1.05)
        mean_at_goal = profile['mean'][bearings.size // 2]
        polar.plot(0, zero_radius + mean_at_goal, 'o', color=color, ms=2.8)
        page_text(center, 338, '0°', ha='center', fontsize=8)
        page_text(center, 395, label, ha='center', color=color, fontsize=8.5)

    d = axes_at(486, 407, 134, plot_bottom-407)
    d.axvline(0, color='#888888', lw=.65, ls=(0,(2,2)), zorder=0)
    for y, level in enumerate(levels[1:]):
        mean, lo, hi = data['full_contrasts'][level]
        d.errorbar(mean, y, xerr=[[mean-lo],[hi-mean]], color=colors[level],
                   fmt='o', ms=4.1, elinewidth=1, capsize=2.2, capthick=1)
    d.set(xlim=(-.34,.04), ylim=(2.57,-.57), xticks=(-.3,-.2,-.1,0),
          yticks=range(3), yticklabels=LABELS[1:],
          xlabel='Change vs. control (per radian)')
    d.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
    setup_axes(d)
    d.spines['left'].set_visible(False)
    d.tick_params(axis='y', length=0, pad=5)
    for text_artist, color in zip(d.get_yticklabels(), PALETTE[1:], strict=True):
        text_artist.set_color(color)
    return fig


def save_figure(fig, output, stem):
    """Save a vector PDF and rasterize it for the review preview."""
    pdf = output / f'{stem}.pdf'
    fig.savefig(pdf, metadata={'Title': stem.replace('_', ' '),
                              'CreationDate': None, 'ModDate': None})
    with pymupdf.open(pdf) as document:
        document[0].get_pixmap(matrix=pymupdf.Matrix(2.7, 2.7), alpha=False).save(
            str(output / f'{stem}.png'))
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path,
                        default=source.ROOT/'.cache/goldshtein-bats-figure'/source.DATA_REVISION)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    parser.add_argument('--font-dir', type=Path)
    parser.add_argument('--font-family', default='Arial')
    parser.add_argument('--photo', type=Path)
    args = parser.parse_args(argv)

    from goldshtein_bats_layout_data import load_layout_data

    fonts = source.resolve_fonts(args.font_dir, args.font_family)
    root = source.prepare_inputs(args.data_root, download=args.download)
    data = load_layout_data(root)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    style = {'font.family': args.font_family, 'font.size': 9.5,
             'text.color': INK, 'axes.labelcolor': INK,
             'axes.edgecolor': INK, 'axes.linewidth': .65,
             'xtick.color': INK, 'ytick.color': INK, 'pdf.fonttype': 42}
    with mpl.rc_context(style):
        save_figure(draw_composite(data, photo_path=args.photo), args.output_dir, 'bat_composite')
    summary = {
        'inputs_sha256': source.INPUT_HASHES,
        'data_revision': source.DATA_REVISION,
        'fonts': fonts,
        'n_bats_by_treatment': data['counts'],
        'n_model_steps': len(data['evidence']),
        'elpd': data['scores'],
        'density_at_zero_by_treatment': data['density_at_zero_by_treatment'],
        'density_at_zero_treatment_contrasts': data['full_contrasts'],
        'goal_window_half_width_degrees': data['goal_window_half_width_degrees'],
        'goal_window_rates': data['goal_window_groups'],
        'goal_window_treatment_contrasts': data['goal_window_contrasts'],
        'goal_window_pooled_contrast': data['pooled_goal_window_contrast'],
    }
    (args.output_dir/'bat_layout_review.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
