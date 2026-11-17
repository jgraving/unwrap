"""Rebuild the six-panel beetle figure from observed bearings and a saved fit.

Run ``python notebooks/beetles_figure.py --help`` from the repository root.
No sampling, notebook execution, browser, or manuscript checkout is needed.
"""
from __future__ import annotations

import argparse
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
CAPTION = (
    "Artificial light shifts dung beetle orientation at multiple levels. "
    "(a) Each of ten beetles repeatedly exits the arena, with ten trials per "
    "lighting condition; four light-on bearings are missing. The arena drawing "
    "is schematic. (b) Curves show kernel density estimates of the posterior "
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


def panel_heading(fig, x, y, letter, title, title_x, *, align="center"):
    fig.text(x, y, letter + ".", fontsize=13, weight="bold", va="top", color=INK)
    fig.text(title_x, y - 0.003, title, fontsize=9.5, weight="bold", va="top",
             ha=align, color=INK, linespacing=1.2)


def draw_polar(fig, arrays, df):
    ax = fig.add_axes([0.61, 0.54, 0.32, 0.375], projection="polar")
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
    with mpl.rc_context(style):
        fig = plt.figure(figsize=(7.2, 4.8), dpi=160)
        panel_heading(fig, 0.025, 0.973, "a", "Repeated exit bearings",
                      0.27)
        panel_heading(fig, 0.545, 0.973, "b", "Reorientation", 0.77)
        fig.text(0.402, 0.832, "Exit\nbearing", ha="center", fontsize=8.5, color=INK, linespacing=1.2)
        lamp_xy = draw_polar(fig, arrays, df)
        fig.text(lamp_xy[0], lamp_xy[1] - 0.052, "Floodlight", ha="center",
                 va="top", fontsize=8.5, color=INK)
        handles = [Line2D([], [], color=c, lw=1.7, marker="o", ms=4, label=t)
                   for c, t in ((OFF, "Light off"), (ON, "Light on"))]
        fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.77, 0.918),
                   ncol=2, frameon=False, fontsize=8, handlelength=1.6, columnspacing=1.7)
        fig.add_artist(Line2D([0.025, 0.975], [0.456, 0.456],
                              transform=fig.transFigure, color=RULE, lw=0.7))
        for x, letter, title, ylabel, key in (
            (0.025, "c", "Within-individual", "Mean resultant length", "within"),
            (0.270, "d", "Between-individual", "Mean resultant length", "between"),
            (0.515, "e", "Floodlight attraction", "Exit probability (±45°)", "mass"),
            (0.760, "f", "Treatment effect", "Change in exit\nprobability (±45°)", "difference"),
        ):
            panel_heading(fig, x, 0.431, letter, title, x + 0.040, align="left")
            draw_summary(fig, x, arrays[key], ylabel,
                         difference=(key == "difference"), uniform=(key == "mass"))
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
    # Original arena artwork, as vector paths only. Coordinates are SVG points.
    arena = ET.parse(ASSET).getroot()
    arena.set("x", "32")
    arena.set("y", "59")
    arena.set("width", "224")
    arena.set("height", "103")
    root.append(arena)
    lamp = ET.parse(FLOODLIGHT_ASSET).getroot()
    width, height = 7.2 * 72, 4.8 * 72
    lamp.set("x", str(lamp_xy[0] * width - 13))
    lamp.set("y", str((1 - lamp_xy[1]) * height - 13))
    lamp.set("width", "26")
    lamp.set("height", "26")
    root.append(lamp)
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
