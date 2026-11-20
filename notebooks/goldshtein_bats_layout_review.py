"""Render the bat composite from the pinned tracks and model fits."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
from matplotlib.transforms import Affine2D
import numpy as np
from PIL import Image
import pymupdf

import goldshtein_bats_publication_figure as source


INK = "#222222"
MAP_FACE = "#f5f6f5"
HEATMAP_NORM = mpl.colors.TwoSlopeNorm(vmin=0, vcenter=.5, vmax=1)
PALETTE = ("#0072B2", "#BC7300", "#008664", "#C95118")
LABELS = ("Control", "Eye cover", "+ Magnet", "+ Anosmia")
PAGE_W, PAGE_H = 646, 489
OUTPUT_DIR = source.ROOT / 'assets' / 'goldshtein-bats-layout-review'


def _directedness_colormap(hues=(286, 129), *, lightness=.71, chroma=.15, neutral_lightness=0):
    """Build two matched, equal-step Oklab ramps from a neutral midpoint."""
    hue = np.deg2rad(hues)
    endpoints = np.column_stack(([lightness] * 2, chroma * np.cos(hue), chroma * np.sin(hue)))
    position = np.linspace(-1, 1, 257)
    neutral = np.array([neutral_lightness, 0, 0])
    lab = neutral + (endpoints[(position >= 0).astype(int)] - neutral) * np.abs(position[:, None])

    # Inverse Oklab transform (https://bottosson.github.io/posts/oklab/).
    lab_to_lms = np.array([
        [1, .3963377774, .2158037573],
        [1, -.1055613458, -.0638541728],
        [1, -.0894841775, -1.2914855480],
    ])
    lms_to_rgb = np.array([
        [4.0767416621, -3.3077115913, .2309699292],
        [-1.2684380046, 2.6097574011, -.3413193965],
        [-.0041960863, -.7034186147, 1.7076147010],
    ])
    linear = (lab @ lab_to_lms.T) ** 3 @ lms_to_rgb.T
    if np.any((linear < 0) | (linear > 1)):
        raise ValueError("Directedness colors lie outside the sRGB gamut")
    rgb = np.where(linear <= .0031308, 12.92 * linear, 1.055 * linear ** (1 / 2.4) - .055)
    return mpl.colors.ListedColormap(rgb, name="directedness").with_extremes(bad="white")


HEATMAP_CMAP = _directedness_colormap()


def _route_shades(color, n_routes):
    """Alternate dark and light tones within one treatment's color."""
    rgb = np.asarray(mpl.colors.to_rgb(color))
    ramp = mpl.colors.LinearSegmentedColormap.from_list(
        'route_shades', (rgb * .48, rgb, rgb + .35 * (1 - rgb))
    )
    order = np.empty(n_routes, dtype=int)
    order[::2] = np.arange((n_routes + 1) // 2)
    order[1::2] = np.arange(n_routes - 1, (n_routes - 1) // 2, -1)
    return ramp(np.linspace(0, 1, n_routes)[order])


def draw_flight_heatmap(ax, scores, *, cmap=HEATMAP_CMAP):
    """Plot one row per individual with fine white row separators."""
    n_bats, n_steps = scores.shape
    step_km = source.STEP_LENGTH_M / 1000
    far, end = (n_steps - .5) * step_km, -.5 * step_km
    ax.imshow(scores, extent=(far, end, n_bats, 0), aspect='auto',
              interpolation='nearest', cmap=cmap, norm=HEATMAP_NORM)
    ax.hlines(np.arange(1, n_bats), end, far, colors='white', linewidth=.42)
    ax.set_xlim(far, end)


def draw_flight_colorbar(fig, cax, *, cmap=HEATMAP_CMAP, tick_fontsize=8):
    """Label the direction of evidence on either side of the neutral midpoint."""
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=HEATMAP_NORM, cmap=cmap),
                      cax=cax, orientation='horizontal')
    cb.solids.set_rasterized(False)
    cb.solids.set_edgecolor('face')
    cb.outline.set_linewidth(.5)
    cb.set_ticks((0, .5, 1), labels=('0\nUndirected', '0.5\nEqual evidence', '1\nDirected'))
    cb.ax.tick_params(labelsize=tick_fontsize, length=2.25, pad=2)
    cb.ax.set_title('Posterior directedness', fontsize=9, pad=3)


def _draw_treatment_icons(ax, drawings, row_ys, colors, *, active_alpha=.54):
    """Repeat the original control glyphs and fade the manipulated senses."""
    # Control-row paths for vision, magnetoreception, olfaction and echolocation.
    glyphs = [drawings[index] for index in (60, 65, 51, 56)]
    faded_columns = ((), (0,), (0, 1), (0, 1, 2))
    for row_y, color, faded in zip(row_ys, colors, faded_columns, strict=True):
        transform = Affine2D().translate(0, row_ys[0] - row_y) + ax.transData
        for column, glyph in enumerate(glyphs):
            x0, _, x1, _ = glyph['rect']
            center = ((x0 + x1) / 2, source.PAGE_H - row_y)
            muted = column in faded
            ax.add_patch(Circle(center, 12, facecolor=color, edgecolor='none',
                                alpha=.15 if muted else active_alpha, clip_on=False, zorder=10))
            patch = source.draw_drawing(ax, glyph, zorder=11)
            patch.set_transform(transform)
            patch.set_alpha(.30 if muted else 1)


def draw_polar_density(ax, bearings, profile, color, radius_max):
    """Draw the full mixture around a common zero-density reference ring."""
    zero_radius = 1.0
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set(ylim=(0, radius_max), xticks=[], yticks=[])
    ax.spines['polar'].set_visible(False)
    ax.plot(bearings, np.full_like(bearings, zero_radius), color='#aaaaaa', lw=.55)
    ax.plot([0, 0], [zero_radius, radius_max], color='#9a9a9a', lw=.6)
    ax.fill_between(bearings, zero_radius + profile['lo'], zero_radius + profile['hi'],
                    color=color, alpha=.26, linewidth=0)
    ax.plot(bearings, zero_radius + profile['mean'], color=color, lw=1.05)
    ax.plot(0, zero_radius + profile['mean'][bearings.size // 2], 'o', color=color, ms=2.8)


def draw_density_contrast(ax, data, colors, *, show_control=True, vertical=False):
    """Plot within-draw changes in mixture density at the goal bearing."""
    levels = data['levels']
    contrast_labels = LABELS if show_control else LABELS[1:]
    contrast_levels = levels if show_control else levels[1:]
    reference_line = ax.axhline if vertical else ax.axvline
    reference_line(0, color='#888888', lw=.65, ls=(0, (2, 2)), zorder=0)
    if show_control:
        ax.plot(0, 0, 'o', color=colors[levels[0]], ms=4.1, zorder=3)
    for index, level in enumerate(levels[1:], start=int(show_control)):
        mean, lo, hi = data['full_contrasts'][level]
        point = (index, mean) if vertical else (mean, index)
        error = {'yerr' if vertical else 'xerr': [[mean-lo], [hi-mean]]}
        ax.errorbar(*point, color=colors[level], fmt='o', ms=4.1, elinewidth=1,
                    capsize=2.2, capthick=1, **error)
    ax.tick_params(length=3, width=.65, pad=3)
    ax.spines[['top', 'right']].set_visible(False)
    if vertical:
        ax.set(xlim=(-.5, len(contrast_labels)-.5), ylim=(-.34, .04),
               xticks=range(len(contrast_labels)), xticklabels=contrast_labels,
               yticks=(-.3, -.2, -.1, 0))
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
        ax.tick_params(axis='x', length=0, pad=5)
        labels = ax.get_xticklabels()
    else:
        ax.set(xlim=(-.34, .04), ylim=(len(contrast_labels)-1+.57, -.57),
               xticks=(-.3, -.2, -.1, 0), yticks=range(len(contrast_labels)),
               yticklabels=contrast_labels)
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0, pad=5)
        labels = ax.get_yticklabels()
    for tick, level in zip(labels, contrast_levels, strict=True):
        tick.set_color(colors[level])


def draw_composite(data, *, photo_path=None, flight_heatmaps=None, heading_scores=None):
    """Compose the photo, matched maps, and three fitted summaries."""
    flight_view = flight_heatmaps is not None
    methods_view = heading_scores is not None
    if methods_view and not flight_view:
        raise ValueError("The methods layout requires flight heatmaps")
    page_h = 565 if flight_view else PAGE_H
    if methods_view:
        page_h = 608
    levels = data["levels"]
    colors = dict(zip(levels, PALETTE, strict=True))
    norm = mpl.colors.Normalize(0, 1)
    cmap = mpl.colormaps["viridis"]
    fig = plt.figure(figsize=(PAGE_W / 72, page_h / 72), facecolor="white")

    def axes_at(x, y, w, h, **kwargs):
        return fig.add_axes([x/PAGE_W, 1-(y+h)/page_h, w/PAGE_W, h/page_h], **kwargs)


    def page_text(x, y, text, **kwargs):
        return fig.text(x/PAGE_W, 1-y/page_h, text, va='center', **kwargs)


    def letter(x, y, value, **kwargs):
        page_text(x, y, value + '.', fontsize=16, weight='bold', **kwargs)


    def setup_axes(ax):
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(length=3, width=.65, pad=3)

    def directedness_scale(x, y, width):
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                          cax=axes_at(x, y, width, 5.5),
                          orientation='horizontal', ticks=[0,.5,1], format='%g')
        cb.solids.set_rasterized(False)
        cb.solids.set_edgecolor('face')
        cb.outline.set_linewidth(.5)
        cb.ax.tick_params(labelsize=8, length=2.25, pad=2)
        page_text(x+width+11, y+2.6, 'Posterior directedness', fontsize=9.5)


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
    icon_rows = (126.58, 158.73, 190.90, 223.05)
    if flight_view:
        _draw_treatment_icons(icon_ax, icons, icon_rows, [colors[level] for level in levels])
    else:
        for index, drawing in sorted(icons.items()):
            patch = source.draw_drawing(icon_ax, drawing, zorder=index+10)
            patch.set_linewidth(patch.get_linewidth()*icon_scale)
    for y, level, label in zip(icon_rows, levels, LABELS, strict=True):
        page_text(108, icon_top+(y-iy0)*icon_scale, label, fontsize=9.3,
                  color=colors[level] if flight_view else INK)

    # b. Four map panels share metric bounds.
    letter(190, 24, 'b')
    map_title = 'Homing routes' if flight_view else 'Homing trajectories'
    page_text(410, 24, map_title, weight='bold', fontsize=11.2, ha='center')
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
    map_h = 128 if methods_view else 216
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
        page_text(center, 45, label, ha='center', weight='bold', fontsize=9.8,
                  color=colors[level] if flight_view else INK)
        for poly in polygons:
            ax.fill(*poly.exterior.xy, fc='#d9e9ed', ec='#a3bcc3', lw=.55, zorder=2)

        unscored, gaps, segments, means = [], [], [], []
        route_colors, gap_colors = [], []
        treatment = tracks.loc[tracks['treatment']==level]
        groups = treatment.groupby('tag', sort=False, observed=True)
        counts[level] = len(groups)
        shades = _route_shades(colors[level], len(groups)) if flight_view else [None] * len(groups)
        for (_, flight), shade in zip(groups, shades, strict=True):
            xy = flight[['x_m','y_m']].to_numpy()
            indexed_segments = np.stack((xy[:-1], xy[1:]), axis=1)
            after_gap = flight['gap_before'].to_numpy()[1:] == 1
            gaps.extend(indexed_segments[after_gap])
            if shade is not None:
                gap_colors.extend([shade] * np.count_nonzero(after_gap))
                route = indexed_segments[~after_gap]
                segments.extend(route)
                route_colors.extend([shade] * len(route))
                continue
            # The heading at a row is the displacement ending at that row.
            score = flight['directedness'].to_numpy()[1:]
            valid = np.isfinite(score) & ~after_gap
            unscored.extend(indexed_segments[~np.isfinite(score) & ~after_gap])
            segments.extend(indexed_segments[valid])
            means.extend(score[valid])

        if flight_view:
            ax.add_collection(LineCollection(segments, colors=route_colors, lw=.70,
                                             alpha=.88, capstyle='round', zorder=3))
        else:
            ax.add_collection(LineCollection(unscored, colors='#a7aeae', lw=.53, alpha=.64, zorder=3))
        ax.add_collection(LineCollection(gaps, colors=gap_colors if flight_view else '#959b9c',
                                         lw=.5, alpha=.65,
                                         linestyles=(0,(2.3,2)), zorder=3))
        if not flight_view:
            colored = LineCollection(segments, cmap=cmap, norm=norm, lw=.76, alpha=.88,
                                     capstyle='round', zorder=4)
            colored.set_array(np.asarray(means))
            ax.add_collection(colored)
        release_size, roost_size = (4.2, 5.8) if methods_view else (5.4, 7.6)
        for site in data['sites']:
            ax.plot(*site, marker='X', ms=release_size, mec='white', mew=.8, mfc=INK, zorder=5)
        ax.plot(*data['roost'], marker='*', ms=roost_size, mec='white', mew=.8, mfc=INK, zorder=5)
        if level == levels[-1]:
            x = xlim[1]-430
            y = ylim[0]+470
            ax.plot([x-1000,x],[y,y], color=INK, lw=1.45, zorder=6)
            ax.annotate('1 km', (x-500,y), textcoords='offset points', xytext=(0,4),
                        ha='center', fontsize=8)
            arrow_y = ylim[1] - (1910 if methods_view else 1470)
            ax.annotate('', (x,arrow_y+700), (x,arrow_y),
                        arrowprops=dict(arrowstyle='-|>', lw=1.05, color=INK,
                                        shrinkA=0, shrinkB=0), zorder=6)
            ax.annotate('N',(x,arrow_y+700), textcoords='offset points', xytext=(0,5),
                        ha='center', fontsize=8)
    assert counts == data['counts'] == dict(zip(levels, (29,14,15,9), strict=True))

    if not flight_view:
        directedness_scale(318, 279, 173)

    # c. Step-level PSIS-LOO ELPD and one standard error per model.
    title_y, plot_top, plot_bottom = (289, 314, 391) if flight_view else (314, 337, 456)
    if methods_view:
        letter(190, 202, 'c')
        page_text(410, 202, 'Model comparison', ha='center', weight='bold', fontsize=10.7)
        page_text(330, 218, 'Absolute ELPD', ha='center', weight='bold', fontsize=9.3)
        comparison = axes_at(270, 232, 120, 50)
    else:
        comparison_y, comparison_top = (309, 330) if flight_view else (title_y, plot_top)
        letter(22, comparison_y, 'c')
        page_text(103, comparison_y, 'Model comparison', ha='center', weight='bold', fontsize=10.7)
        comparison_x, comparison_width = (100, 78) if flight_view else (51, 114)
        comparison = axes_at(comparison_x, comparison_top, comparison_width, plot_bottom-comparison_top)
    for y, key in enumerate(('vm','vmu')):
        score = data['scores'][key]
        comparison.errorbar(score['elpd'], y, xerr=score['se'], color=INK, fmt='o',
                            ms=4.1, elinewidth=1, capsize=2.2, capthick=1, zorder=3)
    lo = min(s['elpd']-s['se'] for s in data['scores'].values())-16
    hi = max(s['elpd']+s['se'] for s in data['scores'].values())+16
    model_names = ('von Mises', 'von Mises\n–Uniform Mixture') if flight_view else ('vM', 'vMU')
    comparison.set(xlim=(lo,hi), ylim=(1.63,-.63), yticks=(0,1),
                   yticklabels=model_names, xlabel='ELPD')
    comparison.xaxis.set_major_locator(MultipleLocator(100))
    setup_axes(comparison)
    comparison.spines['left'].set_visible(False)
    comparison.tick_params(axis='y', length=0, pad=4)
    if flight_view:
        comparison.tick_params(axis='y', labelsize=8.5)
        for label in comparison.get_yticklabels():
            label.set_multialignment('right')

    if methods_view:
        page_text(549, 218, 'Predictive gain by heading error', ha='center',
                  weight='bold', fontsize=9.3)
        gain = axes_at(478, 232, 142, 50)
        gain.axhline(0, color='#888888', linestyle=(0, (2, 2)), lw=.65)
        gain.plot(heading_scores['centers'], heading_scores['means'], 'o', color=INK, ms=3.8)
        gain.set(xlim=(0, 180), xticks=(0, 60, 120, 180),
                 xlabel='Absolute heading error (°)', ylabel='ΔELPD per step')
        gain.yaxis.set_major_locator(MultipleLocator(.1))
        gain.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
        setup_axes(gain)
        gain.tick_params(labelsize=8)
        gain.yaxis.label.set_fontsize(8.5)

    # d. Flight-level or group-level homing progression.
    if flight_view:
        if methods_view:
            title_y, plot_top, plot_bottom = 320, 350, 433
            centers, map_w = (112, 260, 408, 556), 138
            panel_x, panel_center = 22, 334
            x_label_y, colorbar_y = 450, 472
        else:
            panel_x, panel_center = 190, 410
            x_label_y, colorbar_y = 410, 433
        letter(panel_x, title_y, 'd')
        page_text(panel_center, title_y, 'Directedness by individual',
                  ha='center', weight='bold', fontsize=11)
        n_steps = source.LAST_N_STEPS + 1
        row_count = max(counts.values())
        for center, level, treatment in zip(centers, levels, LABELS, strict=True):
            scores = flight_heatmaps[level][:, -n_steps:]
            n_bats = counts[level]
            assert scores.shape == (n_bats, n_steps)
            label = f'{treatment} (n = {n_bats})' if methods_view else f'n = {n_bats}'
            page_text(center, plot_top-7, label, color=colors[level], ha='center', fontsize=8.5)
            ax = axes_at(center-map_w/2, plot_top, map_w, plot_bottom-plot_top)
            draw_flight_heatmap(ax, scores)
            ax.set(ylim=(row_count, 0), yticks=[], xticks=(6, 3, 0))
            ax.spines[['left', 'top', 'right']].set_visible(False)
            ax.tick_params(axis='x', labelsize=8, width=.6, length=2, pad=2)
            if level == levels[0]:
                ax.set_ylabel('Individual', fontsize=9.5, labelpad=4)
        page_text(panel_center, x_label_y, 'Path distance remaining (km)', ha='center', fontsize=9.5)
        draw_flight_colorbar(fig, axes_at(panel_center-66, colorbar_y, 132, 5.5))
    else:
        letter(182, title_y, 'd')
        page_text(309, title_y, 'Homing progression', ha='center', weight='bold', fontsize=11)
        c = axes_at(224, plot_top, 178, plot_bottom-plot_top)
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
    if flight_view:
        offset = 43 if methods_view else 0
        letter_x, density_y = 22, 450 + offset
        polar_centers, polar_top, polar_size = (59, 141, 223, 305), 473 + offset, 72
        contrast_box = (444, 471 + offset, 156, 64)
    else:
        letter_x, density_y = 420, 314
        polar_centers, polar_top, polar_size = (444, 497, 550, 603), 341, 47
        contrast_box = (486, 407, 134, plot_bottom-407)
    letter(letter_x, density_y, 'e')
    page_text(letter_x+23, density_y, 'Predictive heading density', weight='bold', fontsize=11)
    bearings = data['density_bearings']
    profiles = data['density_profiles']
    radius_max = 1.0 + 1.08 * max(profile['hi'].max() for profile in profiles.values())
    for center, level, label in zip(polar_centers, levels, LABELS, strict=True):
        polar = axes_at(center-polar_size/2, polar_top, polar_size, polar_size, projection='polar')
        color = colors[level]
        draw_polar_density(polar, bearings, profiles[level], color, radius_max)
        page_text(center, polar_top-3, '0°', ha='center', fontsize=8)
        page_text(center, polar_top+polar_size+7, label, ha='center', color=color, fontsize=8.5)

    d = axes_at(*contrast_box)
    draw_density_contrast(d, data, colors, show_control=flight_view)
    if flight_view:
        from locust_density_figure import add_compact_unit_label

        x, y, width, height = contrast_box
        prop = d.xaxis.label.get_fontproperties()
        add_compact_unit_label(fig, x+width/2, page_h-(y+height+24),
                               'Density change vs. control', 'radian', -1,
                               fontsize=prop.get_size_in_points(), prop=prop)
    else:
        d.set_xlabel('Change vs. control (per radian)')
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
