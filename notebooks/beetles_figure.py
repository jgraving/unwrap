"""Rebuild the six-panel beetle figure from observed bearings and a saved fit.

Run ``python notebooks/beetles_figure.py --help`` from the repository root.
No sampling, notebook execution, browser, or manuscript checkout is needed.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import platform
import subprocess
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
import xarray as xr

from beetles_plots import (
    OFF, ON, STUDY, STUDY_PATH, FLOOD_DIR, FLOOD_HALFWIDTH,
    _mean_hdi, _vonmises_interval_mass,
)
from beetles_preprocess import load_beetle_model_table
from unwrap.circular import A_of_kappa, circular_mean, resultant_length_angles, von_mises_kde
from unwrap.link.numpy import circular_modulo, softplus

REPO = Path(__file__).resolve().parents[1]
ASSET = Path(__file__).parent / "assets" / "beetle_arena.svg"
FLOODLIGHT_ASSET = ASSET.with_name("beetle_floodlight.svg")
DATA_REVISION = "31eb4bdb4a0fe9a0c0bd20df7f101dbd884a8dcc"
INPUT_HASHES = {
    "beetle_security_light/SecurityLightBergsig.txt":
        "5a3acd579c166c4357204566fbf9eca1233d220af0faf48cffd87ccf2155f78f",
    "fit_cache/beetles_fit.nc":
        "f985f7e55f6a28ffb9d6ec544c624b03d8195b7447ec575a43a82351247bf5d1",
    "fit_cache/beetles_fit.meta.json":
        "51dbcdd857bd26f3e2dadfe90dac85633c2d4698e868cd07330d722f20bbd278",
}
INK = "#262626"
MUTED = "#686868"
RULE = "#dededc"
HALFWIDTH = FLOOD_HALFWIDTH
HDI_PROB = 0.89
KDE_KAPPA = 80.0
SVG = "http://www.w3.org/2000/svg"
ARENA_CENTRE = np.array([151.31826, 77.6316])
ARENA_RADII = np.array([128.976, 50.27067])
ARENA_BOUNDS = {condition: (x, 70, 110, 110 * 122 / 262)
                for condition, x in (("off", 20), ("on", 145))}
SCHEMATIC_PATH_WIDTH = 3.4
CAPTION = (
    "Artificial light shifts dung beetle orientation at multiple levels. "
    "(a) Each of ten beetles repeatedly exits the arena, with ten trials per "
    "lighting condition; four light-on bearings are missing. Matched light-off "
    "and light-on arena schematics show the two conditions experienced by the "
    "same beetles. Illustrative exits are dispersed with the light off and "
    "cluster toward the floodlight with the light on. Paths and beetle "
    "positions are schematic, not observed trajectories or a fitted sample. "
    "(b) Curves show kernel density estimates of the posterior "
    "population mean direction, with a common radial density scale; dots show "
    "each beetle's observed circular mean bearing. Bearings increase clockwise "
    "from north, viewed from above. The lamp icon is schematic. The arc is "
    f"centred at {STUDY['sector_center_degrees']:g}° with a "
    f"{STUDY['sector_halfwidth_degrees']:g}° half-width, using the nominal "
    "south-southwest bearing reported by Foster et al. (2021), rather than a "
    "surveyed azimuth. Sector probabilities and intervals are conditional on "
    "this approximate centre. "
    "(c) Population-level within-"
    "individual concentration, converted to mean resultant length using "
    "I₁(κ)/I₀(κ). (d) Between-individual concentration: the resultant length "
    "of the ten fitted individual mean directions, calculated in each draw. "
    "(e) Probability of an exit within the marked sector, using population "
    "direction and concentration (individual offsets set to zero); the dashed "
    "line marks the uniform expectation of 0.25. (f) The paired draw-wise "
    "difference in this probability, on minus off. Points show posterior "
    "means and bars show 89% highest-density intervals."
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_inputs(root, *, download=False):
    """Fetch only three pinned files, or validate an existing local copy.

    Never replace existing data. A different fit must be selected deliberately
    by updating the revision and checksums together in this source file.
    """
    root = Path(root)
    for relative, expected in INPUT_HASHES.items():
        path = root / relative
        if not path.exists():
            if not download:
                raise FileNotFoundError(f"Missing {path}; supply --download or --data-root.")
            url = f"https://raw.githubusercontent.com/jgraving/unwrap-data/{DATA_REVISION}/{relative}"
            path.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(url, timeout=60) as response:
                content = response.read()
            if hashlib.sha256(content).hexdigest() != expected:
                raise ValueError(f"Downloaded checksum mismatch: {relative}")
            path.write_bytes(content)
        if sha256(path) != expected:
            raise ValueError(f"Checksum mismatch: {path}; expected unwrap-data {DATA_REVISION}.")
    return root


def posterior_arrays(post, *, sector_center=FLOOD_DIR):
    """The estimands used by beetles_plots.plot_concentration, unchanged.

    Use named dimensions so chain/draw storage order cannot change the result.
    This intentionally uses fixed effects for within-individual concentration
    and sector probability, and fitted individual effects for between-individual
    concentration. It does not average predictive densities across beetles.
    """
    names = ("Intercept", "Light", "kappa_Intercept", "kappa_Light")
    for name in (*names, "1|Individual", "Light|Individual"):
        if name not in post:
            raise ValueError(f"Saved posterior is missing {name!r}")
        if not np.isfinite(post[name]).all():
            raise ValueError(f"Non-finite posterior values in {name!r}")
    fixed = {name: post[name].transpose("chain", "draw").values for name in names}
    individual_dim = [d for d in post["1|Individual"].dims if d not in ("chain", "draw")]
    if len(individual_dim) != 1:
        raise ValueError("Expected one individual dimension in the saved posterior")
    dim = individual_dim[0]
    rand_i = post["1|Individual"].transpose("chain", "draw", dim).values
    rand_l = post["Light|Individual"].transpose("chain", "draw", dim).values
    mu = np.stack([
        circular_modulo(fixed["Intercept"]).ravel(),
        circular_modulo(fixed["Intercept"] + fixed["Light"]).ravel(),
    ])
    kappa = np.stack([
        softplus(fixed["kappa_Intercept"]).ravel(),
        softplus(fixed["kappa_Intercept"] + fixed["kappa_Light"]).ravel(),
    ])
    mu_individual = np.stack([
        circular_modulo(fixed["Intercept"][..., None] + rand_i),
        circular_modulo((fixed["Intercept"] + fixed["Light"])[..., None] + rand_i + rand_l),
    ]).reshape(2, -1, rand_i.shape[-1])
    mass = _vonmises_interval_mass(mu, kappa, sector_center - HALFWIDTH, sector_center + HALFWIDTH)
    return {
        "mu": mu,
        "within": A_of_kappa(kappa),
        "between": resultant_length_angles(mu_individual, axis=-1),
        "mass": mass,
        "difference": mass[1] - mass[0],
    }


def resolve_fonts(directory=None, family="Arial"):
    """Require the requested font; never silently substitute another family."""
    if directory:
        fonts = font_manager.findSystemFonts(fontpaths=[str(directory)])
        for font in fonts:
            font_manager.fontManager.addfont(font)
    resolved = {}
    for weight in ("normal", "bold"):
        try:
            path = font_manager.findfont(
                font_manager.FontProperties(family=family, weight=weight),
                fallback_to_default=False,
            )
        except ValueError as exc:
            raise ValueError(
                f"{family} is unavailable. Use --font-dir with licensed font files, "
                "or explicitly choose another family with --font-family."
            ) from exc
        resolved[weight] = {"file": Path(path).name, "sha256": sha256(path)}
    return resolved


def panel_heading(fig, x, y, letter, title, title_x, *, align="center", letter_gap=None):
    heading = fig.text(title_x, y - 0.003, title, fontsize=9.5, weight="bold", va="top",
                       ha=align, color=INK, linespacing=1.2)
    letter_align = "left"
    if letter_gap is not None:
        # Keep the heading centred on its axes and place its letter a fixed
        # number of typographic points before the rendered title.
        box = heading.get_window_extent(renderer=fig.canvas.get_renderer())
        left = fig.transFigure.inverted().transform((box.x0, 0))[0]
        x = left - letter_gap / (72 * fig.get_figwidth())
        letter_align = "right"
    fig.text(x, y, letter + ".", fontsize=13, weight="bold", va="top",
             ha=letter_align, color=INK)


def placed_vector(asset, prefix, *, x, y, width, height, path_color=None):
    """Place the original path-only artwork with unique IDs for each copy.

    The four arena exit paths use a common, legible stroke width.
    """
    root = ET.parse(asset).getroot()
    for element in root.iter():
        original_id = element.get("id")
        if path_color is not None and original_id in {"path3", "path4", "path23", "path24"}:
            element.set("style", element.get("style")
                        .replace("stroke:#3e1f51", f"stroke:{path_color}")
                        .replace("stroke-width:1.33333", f"stroke-width:{SCHEMATIC_PATH_WIDTH}")
                        .replace("stroke-linecap:butt", "stroke-linecap:round")
                        .replace("stroke-linejoin:miter", "stroke-linejoin:round"))
        if original_id:
            element.set("id", f"{prefix}-{original_id}")
    root.set("id", prefix)
    for key, value in (("x", x), ("y", y), ("width", width), ("height", height)):
        root.set(key, str(value))
    root.set("overflow", "visible")
    return root


def svg_matrix(linear, translation):
    values = (*linear[:, 0], *linear[:, 1], *translation)
    return "matrix(" + " ".join(f"{v:.9g}" for v in values) + ")"


def rolling_beetle(original, prefix, index, position, radial, size):
    """Reuse four original views with raised hind legs resting on each ball.

    Keep the camera's upright direction: choose/mirror the closest view and
    limit its turn to 25 degrees. Place the body and ball separately, drawing
    the ball in front and two filled, tapered hind legs over its upper cap.
    """
    poses = (
        dict(ball=14, body=(22, 21, 20), ground=(15, 16, 17),
             anchor=(221.41477, 115.67518), head=(199.5, 102), flip=(1, 1), lift=1,
             hips=((203.35827, 104.24067), (213.2612, 101.5208))),
        dict(ball=25, body=(26, 32, 33), ground=(27, 28, 29),
             anchor=(151.07448, 22.453), head=(150, 36), flip=(1, -1), lift=8,
             hips=((145.1332, 30.884533), (156.314, 31.298667))),
        dict(ball=34, body=(35, 41, 42), ground=(37, 38, 39),
             anchor=(82.74406, 32.1399), head=(87.5, 45), flip=(1, -1), lift=12,
             hips=((82.757467, 40.500533), (91.170133, 38.514533))),
        dict(ball=5, body=(8, 7, 6), ground=(9, 10, 11),
             anchor=(239.24303, 38.07613), head=(210.5, 40.5), flip=(-1, 1), lift=0,
             hips=((226.02587, 31.765467), (223.12053, 28.860133))),
    )
    pose = poses[index]
    mirror = np.diag(pose["flip"])
    source = mirror @ (np.array(pose["anchor"]) - np.array(pose["head"]))
    angle = np.angle(np.exp(1j * (np.arctan2(radial[1], radial[0])
                                  - np.arctan2(source[1], source[0]))))
    angle = np.clip(angle, -np.deg2rad(25), np.deg2rad(25))
    linear = size * np.array([[np.cos(angle), -np.sin(angle)],
                              [np.sin(angle), np.cos(angle)]]) @ mirror
    ball_shift = position - linear @ np.array(pose["anchor"])
    body_shift = ball_shift - np.array([0, pose["lift"] * size])
    beetle = ET.Element(f"{{{SVG}}}g", id=f"{prefix}-beetle-{index}",
                        attrib={"data-source-pose": str(pose["ball"])})

    def add_part(number, shift):
        part = deepcopy(original[f"{prefix}-path{number}"])
        part.set("id", f"{prefix}-beetle-{index}-path{number}")
        part.set("transform", svg_matrix(linear, shift))
        beetle.append(part)

    for number in (*pose["ground"], *pose["body"]):
        add_part(number, body_shift)
    add_part(pose["ball"], ball_shift)
    hips = sorted((linear @ np.array(hip) + body_shift for hip in pose["hips"]),
                  key=lambda point: point[0])
    for side, hip in zip((-1, 1), hips):
        contact = position + size * np.array([side * 5.0, -8.5])
        knee = (hip + contact) / 2 + size * np.array([side * 2.8, -3.5])
        toe = contact + size * np.array([side * 0.7, 1.3])
        centres = np.array([hip, knee, contact, toe])
        segments = np.diff(centres, axis=0)
        normals = np.column_stack((-segments[:, 1], segments[:, 0]))
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        joints = np.vstack((normals[0], normals[:-1] + normals[1:], normals[-1]))
        joints /= np.linalg.norm(joints, axis=1)[:, None]
        # Broad articulated upper leg, narrowing through the shank to the
        # foot. Mitred edges match the angular filled legs in the source SVG.
        half_widths = size * np.array([0.95, 1.15, 0.42, 0.08])
        mitres = np.ones(4)
        mitres[1:-1] = np.clip(np.sum(joints[1:-1] * normals[:-1], axis=1), 0.5, 1)
        offsets = joints * (half_widths / mitres)[:, None]
        outline = np.vstack((centres + offsets, (centres - offsets)[::-1]))
        points = [" ".join(f"{value:.8g}" for value in point) for point in outline]
        ET.SubElement(beetle, f"{{{SVG}}}path", id=f"{prefix}-beetle-{index}-hind-{side}",
                      d="M " + " L ".join(points) + " Z",
                      style=f"fill:#575756;stroke:#3b493e;stroke-width:{0.548 * size:.8g};"
                            "stroke-linecap:butt;stroke-linejoin:miter")
    return beetle, linear @ np.array(pose["head"]) + body_shift


def arena_vector(condition):
    """Compose the treatment schematic from original vector paths.

    The on exits fan toward the nominal lamp bearing on the arena ellipse.
    Four original beetle views are adapted separately, with rear legs on the
    upper ball surfaces. Size and drawing order vary with foreground depth.
    """
    prefix = f"arena-light-{condition}"
    x, y, width, height = ARENA_BOUNDS[condition]
    root = placed_vector(ASSET, prefix, x=x, y=y, width=width, height=height,
                         path_color=OFF if condition == "off" else ON)
    if condition == "off":
        return root
    original = {element.get("id"): element for element in root}
    ground = original[f"{prefix}-path2"]
    title = root.find(f"{{{SVG}}}title")
    root[:] = [title, ground]
    # Original path IDs and foreground sprite anchors in SVG coordinates.
    path_ids = ("path3", "path4", "path23", "path24")
    bearings = FLOOD_DIR + np.deg2rad([-40, -12, 13, 40])
    beetles = []
    for index, (path_id, bearing) in enumerate(zip(path_ids, bearings)):
        rim = ARENA_CENTRE + ARENA_RADII * np.array([np.sin(bearing), -np.cos(bearing)])
        radial = rim - ARENA_CENTRE
        direction = radial / np.linalg.norm(radial)
        depth = (rim[1] - ARENA_CENTRE[1]) / ARENA_RADII[1]
        size = 0.94 + 0.08 * depth
        beetle, trail_end = rolling_beetle(original, prefix, index, rim, radial, size)
        # Two gentle cubic curves reach the back of each beetle. Monotone
        # progress along the exit ray avoids loops when the path is short.
        target = trail_end - ARENA_CENTRE
        normal = np.array([-direction[1], direction[0]])
        bend = 3.2 if index % 2 == 0 else -3.2
        points = [ARENA_CENTRE + fraction * target + offset * normal
                  for fraction, offset in ((0, 0), (0.15, bend), (0.35, -bend),
                                           (0.5, 0), (0.65, bend), (0.85, -bend), (1, 0))]
        coordinates = [" ".join(f"{value:.8g}" for value in point) for point in points]
        path = deepcopy(original[f"{prefix}-{path_id}"])
        path.set("d", f"M {coordinates[0]} C {' '.join(coordinates[1:4])} C {' '.join(coordinates[4:])}")
        root.append(path)
        beetles.append((rim[1], beetle))
    for _, beetle in sorted(beetles, key=lambda item: item[0]):
        root.append(beetle)
    return root


def schematic_lamp_position():
    """Project the nominal lamp bearing onto the same ellipse as the exits."""
    point = ARENA_CENTRE + 2.05 * ARENA_RADII * np.array([np.sin(FLOOD_DIR), -np.cos(FLOOD_DIR)])
    x, y, width, _ = ARENA_BOUNDS["on"]
    return np.array([x, y]) + (point - np.array([20, 8])) * width / 262


def draw_polar(fig, arrays, df):
    ax = fig.add_axes([0.66, 0.54, 0.32, 0.375], projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_axis_off()
    theta = np.linspace(-np.pi, np.pi, 721)
    density = np.stack([von_mises_kde(samples, theta, KDE_KAPPA) for samples in arrays["mu"]])
    scale = density.max()
    base, amplitude, ring = 0.40, 0.48, 1.0
    for angle in np.arange(4) * np.pi / 2:
        ax.plot([angle, angle], [0, ring], lw=0.5, color=RULE, zorder=0)
    ax.plot(theta, np.full_like(theta, ring), lw=0.7, color="#c9c9c9", zorder=1)
    for condition, color in enumerate((OFF, ON)):
        radius = base + amplitude * density[condition] / scale
        ax.fill_between(theta, base, radius, color=color, alpha=0.17, linewidth=0)
        ax.plot(theta, radius, color=color, lw=1.7, solid_capstyle="round")
        for _, group in df[df.Light == condition].groupby("Individual", observed=True):
            bearing = circular_mean(group.HeadingRadians.to_numpy())
            ax.scatter([bearing], [ring], s=26, color=color, edgecolor="white", linewidth=0.45, zorder=5)
    # Clockwise compass bearings, viewed from above. Move S outside the
    # marked sector; the other labels sit just beyond the observation ring.
    for angle, label in zip(np.arange(4) * np.pi / 2, ("N", "E", "S", "W")):
        radius = 1.31 if label == "S" else 1.13
        ax.text(angle, radius, label, ha="center", va="center",
                color=MUTED, fontsize=7.5)
    sector = np.linspace(FLOOD_DIR - HALFWIDTH, FLOOD_DIR + HALFWIDTH, 100)
    ax.plot(sector, np.full_like(sector, 1.15), color=INK, lw=0.9)
    for angle in (sector[0], FLOOD_DIR, sector[-1]):
        ax.plot([angle, angle], [1.12, 1.18], color=INK, lw=0.9)
    ax.set_ylim(0, 1.36)
    # Polar axes shrink their box to remain circular at draw time. Apply that
    # aspect adjustment before converting the central spoke into figure
    # coordinates, so the lamp stays aligned with the rendered sector tick.
    ax.apply_aspect()
    point = ax.transData.transform((FLOOD_DIR, 1.39))
    return fig.transFigure.inverted().transform(point)


def draw_summary(fig, x, samples, ylabel, *, difference=False, uniform=False):
    ax = fig.add_axes([x + 0.079, 0.105, 0.135, 0.27])
    ax.set_xlim(-0.52, 1.52)
    lower = min(0.0, min(_mean_hdi(samples, HDI_PROB)) - 0.025) if difference else 0.0
    ax.set_ylim(lower, 1.025)
    ticks = [0, 0.25, 0.5, 0.75, 1]
    if lower < -0.25:
        ticks = list(np.arange(-0.25, lower, -0.25)[::-1]) + ticks
    ax.set_yticks(ticks, [f"{v:g}" for v in ticks])
    ax.tick_params(axis="both", length=3, width=0.7, direction="out", pad=4,
                   labelsize=8, colors=INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_visible(True)
        ax.spines[side].set_linewidth(0.7)
        ax.spines[side].set_color(INK)
    ax.spines["left"].set_bounds(lower, 1)
    ax.spines["bottom"].set_position(("data", 0))
    ax.grid(False)
    ax.set_ylabel(ylabel, fontsize=8, color=INK, labelpad=7, linespacing=1.25)
    if uniform:
        ax.axhline(0.25, color=MUTED, lw=0.9, linestyle=(0, (3, 3)), zorder=1)
        ax.text(1.49, 0.22, "Uniform", ha="right", va="top", color=INK, fontsize=7.5)
    series = [(0.5, samples, ON)] if difference else [(0, samples[0], OFF), (1, samples[1], ON)]
    for location, draws, color in series:
        mean, lo, hi = _mean_hdi(draws, hdi_prob=HDI_PROB)
        # Draw the interval directly: a mean need not lie inside a shortest HDI.
        ax.plot([location, location], [lo, hi], color=color, lw=1.5, zorder=3)
        ax.plot([location - 0.075, location + 0.075], [lo, lo], color=color, lw=1.5)
        ax.plot([location - 0.075, location + 0.075], [hi, hi], color=color, lw=1.5)
        ax.scatter([location], [mean], s=32, color=color, edgecolor="white", linewidth=0.45, zorder=4)
    ax.set_xticks([0.5] if difference else [0, 1], ["On − off"] if difference else ["Off", "On"])
    return ax


def build_svg(arrays, df, *, family="Arial"):
    """Compose a vector figure; text outlines preserve the exact chosen font."""
    style = {
        "font.family": family, "font.size": 8, "text.color": INK,
        "svg.fonttype": "path", "svg.hashsalt": "unwrap-beetles-publication-v1",
        "figure.facecolor": "white", "axes.facecolor": "white",
    }
    arena_lamp_xy = schematic_lamp_position()
    with mpl.rc_context(style):
        fig = plt.figure(figsize=(7.2, 4.8), dpi=160)
        panel_heading(fig, 0.025, 0.973, "a", "Repeated exit bearings",
                      0.27)
        for centre, label in ((75 / 518.4, "Light off"), (200 / 518.4, "Light on")):
            fig.text(centre, 0.845, label, ha="center", va="top",
                     fontsize=9, weight="bold", color=INK)
        fig.text(arena_lamp_xy[0] / 518.4, 1 - (arena_lamp_xy[1] + 16) / 345.6,
                 "Floodlight", ha="center",
                 va="top", fontsize=8, color=INK)
        panel_heading(fig, 0.545, 0.973, "b", "Reorientation", 0.82)
        lamp_xy = draw_polar(fig, arrays, df)
        fig.text(lamp_xy[0], lamp_xy[1] - 0.052, "Floodlight", ha="center",
                 va="top", fontsize=8.5, color=INK)
        handles = [Line2D([], [], color=c, marker="s", linestyle="none", ms=6,
                          markeredgewidth=0, label=t)
                   for c, t in ((OFF, "Light off"), (ON, "Light on"))]
        fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.82, 0.918),
                   ncol=2, frameon=False, fontsize=8, handlelength=0.9,
                   handletextpad=0.5, columnspacing=1.7)
        levels = [
            Line2D([], [], color=INK, lw=1.7, label="Population-level\nmean"),
            Line2D([], [], color=INK, lw=0, marker="o", ms=4,
                   label="Individual-level\nmeans"),
        ]
        fig.legend(handles=levels, loc="upper left", bbox_to_anchor=(0.535, 0.82),
                   frameon=False, fontsize=8, handlelength=1.6, labelspacing=1.4,
                   borderaxespad=0, handletextpad=0.8)
        for x, letter, title, ylabel, key in (
            (0.025, "c", "Within-individual", "Mean resultant length", "within"),
            (0.270, "d", "Between-individual", "Mean resultant length", "between"),
            (0.515, "e", "Floodlight attraction", "Exit probability (±45°)", "mass"),
            (0.760, "f", "Treatment effect", "Change in exit\nprobability (±45°)", "difference"),
        ):
            ax = draw_summary(fig, x, arrays[key], ylabel,
                              difference=(key == "difference"), uniform=(key == "mass"))
            box = ax.get_position()
            panel_heading(fig, x, 0.431, letter, title, box.x0 + box.width / 2,
                          letter_gap=4)
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg", metadata={"Date": None, "Creator": "unwrap beetles_figure.py"})
        plt.close(fig)
    ET.register_namespace("", SVG)
    root = ET.fromstring(buffer.getvalue())
    root.set("role", "img")
    root.set("aria-labelledby", "figure-title figure-description")
    title = ET.Element(f"{{{SVG}}}title", id="figure-title")
    title.text = "Artificial light shifts dung beetle orientation at multiple levels"
    desc = ET.Element(f"{{{SVG}}}desc", id="figure-description")
    desc.text = CAPTION
    root.insert(0, title)
    root.insert(1, desc)
    # Matched treatment schematics from the same original vector artwork.
    # Coordinates are SVG points; paths illustrate the treatment result.
    for condition in ("off", "on"):
        root.append(arena_vector(condition))
    root.append(placed_vector(FLOODLIGHT_ASSET, "arena-floodlight",
                              x=arena_lamp_xy[0] - 10, y=arena_lamp_xy[1] - 10,
                              width=20, height=20))
    width, height = 7.2 * 72, 4.8 * 72
    root.append(placed_vector(FLOODLIGHT_ASSET, "polar-floodlight",
                              x=lamp_xy[0] * width - 13,
                              y=(1 - lamp_xy[1]) * height - 13, width=26, height=26))
    return "\n".join(line.rstrip() for line in ET.tostring(root, encoding="unicode").splitlines()) + "\n"


def write_outputs(output, arrays, df, metadata, *, family="Arial"):
    import cairosvg

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    stem = output / "beetles_publication"
    svg = build_svg(arrays, df, family=family)
    stem.with_suffix(".svg").write_text(svg, encoding="utf-8")
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(stem.with_suffix(".pdf")))
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(stem.with_suffix(".png")), output_width=1800)
    # A self-contained review page; the scientific figure remains pure vector.
    html = f'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Beetle orientation figure</title>
<style>body{{margin:0;background:#f5f5f4;color:#262626;font:15px/1.65 Arial,sans-serif}}
main{{max-width:1100px;margin:40px auto;padding:0 24px}}figure{{margin:0}}
.figure{{background:white;padding:20px;border:1px solid #e6e6e6}}svg{{display:block;width:100%;height:auto}}
figcaption{{max-width:94ch;margin:24px auto}}a{{color:#665017}}nav{{margin-bottom:18px}}</style>
<main><nav><a href="beetles_publication.svg">SVG</a> · <a href="beetles_publication.pdf">PDF</a></nav>
<figure><div class="figure">{svg}</div><figcaption>{CAPTION}</figcaption></figure></main></html>'''
    stem.with_suffix(".html").write_text(html, encoding="utf-8")
    summary = {}
    for key in ("within", "between", "mass", "difference"):
        draws = arrays[key][None, :] if key == "difference" else arrays[key]
        labels = ("on_minus_off",) if key == "difference" else ("off", "on")
        summary[key] = {
            label: dict(zip(("mean", "hdi_low", "hdi_high"), _mean_hdi(values, HDI_PROB)))
            for label, values in zip(labels, draws)
        }
    metadata["summaries"] = summary
    metadata["outputs_sha256"] = {suffix: sha256(stem.with_suffix(suffix)) for suffix in (".svg", ".pdf", ".png", ".html")}
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return stem


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=REPO / ".cache" / "beetle-figure" / DATA_REVISION)
    parser.add_argument("--download", action="store_true", help="Fetch missing pinned beetle inputs (about 4 MB).")
    parser.add_argument("--output-dir", type=Path, default=REPO / "figures" / "beetles")
    parser.add_argument("--font-dir", type=Path, help="Directory containing licensed Arial regular/bold fonts.")
    parser.add_argument("--font-family", default="Arial", help="Explicit font override; default: Arial.")
    args = parser.parse_args(argv)
    fonts = resolve_fonts(args.font_dir, args.font_family)
    data_root = prepare_inputs(args.data_root, download=args.download)
    df = load_beetle_model_table(data_root=data_root)
    fit_path = data_root / "fit_cache" / "beetles_fit.nc"
    with xr.open_dataset(fit_path, group="posterior", engine="h5netcdf") as dataset:
        post = dataset.load()
    with xr.open_dataset(fit_path, group="sample_stats", engine="h5netcdf") as stats:
        divergences = int(stats.diverging.sum())
    if divergences:
        raise ValueError(f"Saved fit has {divergences} divergent transitions; review before plotting.")
    arrays = posterior_arrays(post)
    observed_ids = {str(value) for value in df.Individual.unique()}
    fitted_ids = set(post["Individual__factor_dim"].values.astype(str))
    if observed_ids != fitted_ids or df.groupby("Light").size().to_dict() != {0: 100, 1: 96}:
        raise ValueError("The pinned dataset does not match the ten-beetle figure design.")
    try:
        code_revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        code_revision = None
    metadata = {
        "data_repository": "jgraving/unwrap-data", "data_revision": DATA_REVISION,
        "inputs_sha256": INPUT_HASHES, "code_revision": code_revision,
        "source_sha256": {p.relative_to(REPO).as_posix(): sha256(p) for p in
                          (Path(__file__), ASSET, FLOODLIGHT_ASSET, STUDY_PATH, Path(__file__).with_name("beetles_plots.py"),
                           Path(__file__).with_name("beetles_preprocess.py"),
                           REPO / "unwrap" / "circular.py", REPO / "unwrap" / "bessel.py",
                           REPO / "unwrap" / "link" / "numpy.py")},
        "python": platform.python_version(),
        "packages": {name: importlib.metadata.version(name) for name in
                     ("numpy", "scipy", "pandas", "matplotlib", "arviz", "xarray", "h5netcdf", "h5py", "cairosvg")},
        "font_family": args.font_family, "fonts": fonts,
        "posterior_chains": post.sizes["chain"], "posterior_draws_per_chain": post.sizes["draw"],
        "divergences": divergences, "observations": {"off": 100, "on": 96},
        "hdi_probability": HDI_PROB, "kde_kappa": KDE_KAPPA,
        **STUDY,
        "caption": CAPTION,
    }
    stem = write_outputs(args.output_dir, arrays, df, metadata, family=args.font_family)
    print(f"Wrote SVG, PDF, PNG, HTML and provenance JSON: {stem}")


if __name__ == "__main__":
    main()
