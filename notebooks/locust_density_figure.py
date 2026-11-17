"""Rebuild the locust density/order publication figure without sampling.

Run from the repository root; see LOCUST_DENSITY_FIGURE.md for provenance,
the corrected predictor coordinates, and source-image extraction.
"""
from __future__ import annotations

import argparse
import base64
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
from matplotlib.patches import Rectangle
import numpy as np
from scipy.special import i0e
import xarray as xr

from unwrap.circular import circular_mean
from unwrap.bessel import inverse_bessel_ratio
from unwrap.link.numpy import circular_modulo, softplus

ROOT = Path(__file__).resolve().parents[1]
ASSETS = Path(__file__).parent / "assets/locust_density_order"
DATA_REVISION = "31eb4bdb4a0fe9a0c0bd20df7f101dbd884a8dcc"
INPUT_HASHES = {
    "fit_cache/locust_density_hier.nc": "af8df700a42354c59fa673e8635d25e4c11f038f404cdb9513aaf5b853176498",
    "fit_cache/locust_density_hier.meta.json": "ebbbd2a67f5145961da05357a32c1cb9f48181b37bd598050a48c4efaeaeec49",
}
INK, MUTED = "#262626", "#686868"
MU_COLOR, KAPPA_COLOR = "#8c564b", "#2ca02c"
ORDER_LEVELS = np.array([0, .2, .4, .6, .8, 1.])
DENSITY_LEVELS = np.array([1, 4, 8, 16, 32, 64.])
TERMS = ("Intercept", "t_density", "t_order", "t_density:t_order")
SVG = "http://www.w3.org/2000/svg"
CAPTION = (
    "Order and density of a virtual swarm shape locust directional responses. "
    "(a) Real locust. (b) Original front, lateral, and rear renderings of the "
    "virtual locust. (c) Panoramic VR arena and spherical treadmill, with the "
    "projection screen raised. (d) Real and virtual swarms; arrow indicates "
    "the virtual march direction. (e) Seeded illustrations of density and "
    "motion order; these are stimulus schematics, not recorded trajectories. "
    "(f) Five screening specifications fitted to 440 recording mean directions; "
    "D denotes density and O denotes order. "
    "points show PSIS-LOO ELPD differences from the best model, bars ±1 standard "
    "error. Values reproduce the versioned Python notebook's saved screening "
    "table, printed to six decimals. (g) Population fixed-effect coefficients "
    "from the saved hierarchical fit to 18,641 headings from 440 recordings. "
    "Points show posterior means (circular for the direction intercept), bars "
    "95% highest-density intervals. Density is log₂-scaled from 1–64 m⁻² to "
    "[−1, 1]; order is scaled from 0–1 to [−1, 1]. Coefficients are on the "
    "pre-link scale. (h) Reference density at the march direction, f(0), with "
    "the recording concentration offset set to zero. Curves are exp(mean(log "
    "f(0))), with exponentiated 90% highest-density intervals of log f(0). "
    "These are geometric posterior summaries of density per radian, not "
    "probabilities or ELPD scores. Dashed lines mark uniform density 1/(2π). "
    "Curves use the original scaled coordinates of the stored Bambi "
    "coefficients, correcting the notebook helper's extra centering. "
    "Photos and stimulus renderings are extracted from the manuscript figure; "
    "photo credit: Centre for the Advanced Study of Collective Behaviour "
    "(manuscript states used with permission)."
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_inputs(root, *, download=False):
    root = Path(root)
    for relative, expected in INPUT_HASHES.items():
        path = root / relative
        if not path.exists():
            if not download:
                raise FileNotFoundError(f"Missing {path}; supply --download or --data-root")
            url = f"https://raw.githubusercontent.com/jgraving/unwrap-data/{DATA_REVISION}/{relative}"
            with urllib.request.urlopen(url, timeout=60) as response:
                data = response.read()
            if hashlib.sha256(data).hexdigest() != expected:
                raise ValueError(f"Downloaded checksum mismatch: {relative}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        if sha256(path) != expected:
            raise ValueError(f"Checksum mismatch: {path}; existing files are never overwritten")
    return root


def resolve_fonts(directory=None, family="Arial"):
    if directory:
        for path in font_manager.findSystemFonts(fontpaths=[str(directory)]):
            font_manager.fontManager.addfont(path)
    fonts = {}
    for weight in ("normal", "bold"):
        try:
            path = font_manager.findfont(font_manager.FontProperties(family=family, weight=weight),
                                         fallback_to_default=False)
        except ValueError as exc:
            raise ValueError(f"{family} unavailable; supply licensed fonts with --font-dir") from exc
        fonts[weight] = {"file": Path(path).name, "sha256": sha256(path)}
    return fonts


def hdi(samples, probability):
    """Shortest contiguous interval, matching the notebook's ArviZ convention."""
    values = np.sort(np.asarray(samples), axis=0)
    count = int(np.floor(probability * len(values)))
    widths = values[count:] - values[:len(values) - count]
    indices = np.argmin(widths, axis=0)
    return np.take_along_axis(values, indices[None], axis=0)[0], np.take_along_axis(values, (indices + count)[None], axis=0)[0]


def reference_parameters(post, order, density):
    """Stored Bambi intercepts have already been uncentred after sampling.

    Use the original scaled predictors and their product. In particular, do
    not subtract training means again or multiply centred marginal columns.
    The concentration random intercept is zero for these reference curves.
    """
    order, density = np.broadcast_arrays(np.asarray(order, float), np.asarray(density, float))
    if (not np.isfinite(order).all() or not np.isfinite(density).all()
            or np.any((order < 0) | (order > 1) | (density < 1) | (density > 64))):
        raise ValueError("Prediction requested outside the fitted stimulus range")
    o, d = 2 * order.ravel() - 1, np.log2(density.ravel()) / 3 - 1
    design = np.stack([np.ones_like(o), d, o, d * o])
    eta = []
    for prefix in ("", "kappa_"):
        coefs = np.stack([post[prefix + term].transpose("chain", "draw").values.ravel() for term in TERMS])
        if not np.isfinite(coefs).all():
            raise ValueError("Non-finite posterior coefficients")
        eta.append(coefs.T @ design)
    return circular_modulo(eta[0]), softplus(eta[1])


def reference_summary(post, order, density):
    mu, kappa = reference_parameters(post, order, density)
    # Stable log density at zero: i0(k) = exp(k) * i0e(k).
    logp = kappa * (np.cos(mu) - 1) - np.log(2 * np.pi) - np.log(i0e(kappa))
    low, high = hdi(logp, .90)
    return {"mean": np.exp(logp.mean(axis=0)), "low": np.exp(low), "high": np.exp(high)}


def statistics(post):
    order_grid = np.linspace(0, 1, 101)
    density_grid = np.geomspace(1, 64, 101)
    coefficients = []
    for term in TERMS:
        row = []
        for prefix in ("", "kappa_"):
            values = post[prefix + term].values.ravel()
            if term == "Intercept" and prefix == "":
                mean = float(circular_mean(values))
                lo, hi = hdi(circular_modulo(values - mean), .95)
                lo, hi = lo + mean, hi + mean
            else:
                mean = float(values.mean())
                lo, hi = hdi(values, .95)
            row.append([mean, float(lo), float(hi)])
        coefficients.append(row)
    return {"order_grid": order_grid, "density_grid": density_grid,
            "over_order": [reference_summary(post, order_grid, d) for d in DENSITY_LEVELS],
            "over_density": [reference_summary(post, o, density_grid) for o in ORDER_LEVELS],
            "coefficients": np.asarray(coefficients)}


def style_axis(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(.65)
        ax.spines[side].set_color(INK)
    ax.tick_params(length=3, width=.65, pad=3, colors=INK, labelsize=7)
    ax.grid(False)


def heading(fig, letter, title, x, y, center):
    title_artist = fig.text(center, y - .002, title, ha="center", va="top",
                            fontsize=9.2, weight="bold")
    label = fig.text(x, y, letter + ".", va="top", fontsize=12, weight="bold")
    renderer = fig.canvas.get_renderer()
    title_box = title_artist.get_window_extent(renderer)
    label_box = label.get_window_extent(renderer)
    label.set_x((title_box.x0 - label_box.width - 4 * fig.dpi / 72) / fig.bbox.width)


def photo(fig, name, bounds, *, crop=None, fit="cover", rotation=0):
    # Place native JPEGs into the final SVG directly, avoiding Matplotlib's
    # image resampling. A transparent axes placeholder records the placement.
    ax = fig.add_axes(bounds)
    ax.set_axis_off()
    return dict(name=name, bounds=bounds, crop=crop, fit=fit, rotation=rotation)


def stimulus_schematic(fig):
    left, bottom, width, height = .281, .125, .174, .484
    ax = fig.add_axes([left, bottom, width, height])
    ax.set_axis_off(); ax.set_xlim(-.35, 2.6); ax.set_ylim(-.2, 6.2)
    rng = np.random.default_rng(123)
    for i, (density, order) in enumerate(zip(DENSITY_LEVELS[::-1], ORDER_LEVELS[::-1])):
        y = 5 - i
        red, purple = plt.cm.Reds(.95 - .10 * i), plt.cm.Purples(.95 - .10 * i)
        boxes = []
        for x in (0, 1.3):
            box = Rectangle((x, y), 1, .82, facecolor="#fafafa", edgecolor="#dededc", lw=.6)
            ax.add_patch(box); boxes.append(box)
        xy = rng.uniform([.07, y + .07], [.93, y + .75], (int(density), 2))
        ax.scatter(xy[:, 0], xy[:, 1], s=6, color=red, linewidths=0)
        for start in np.linspace(.1, .72, 7):
            kappa = inverse_bessel_ratio(float(min(order, .999)))
            angles = rng.vonmises(0, kappa, 120)
            increments = np.column_stack([np.cos(angles), np.sin(angles)])
            path = np.vstack([[0, 0], np.cumsum(increments, axis=0)])
            span = np.maximum(np.ptp(path, axis=0), 1e-9)
            scale = min(.92 / span[0], .68 / span[1])
            path = (path - (path.max(axis=0) + path.min(axis=0))/2) * scale
            offset = (start - .41) * max(0, 1 - np.ptp(path[:, 1])/.68)
            line, = ax.plot(1.8 + path[:, 0], y + .41 + offset + path[:, 1],
                            color=purple, lw=.8)
            line.set_clip_path(boxes[1])
        ax.text(-.12, y + .41, f"{density:g}", ha="right", va="center", fontsize=7.5)
        ax.text(2.42, y + .41, f"{order:g}", va="center", fontsize=7.5)
    ax.text(.5, 6.45, "Density\n(m$^{-2}$)", ha="center", va="top", fontsize=7.5)
    ax.text(1.8, 6.45, "Order", ha="center", va="top", fontsize=7.5)
    ax.annotate("", xy=(2.2, 6.12), xytext=(1.4, 6.12),
                arrowprops={"arrowstyle": "-|>", "lw": .7, "color": INK})


def draw_figure(stats, selection, *, family="Arial"):
    with mpl.rc_context({"font.family": family, "font.size": 7.5, "text.color": INK,
                         "mathtext.fontset": "custom", "mathtext.rm": family,
                         "mathtext.it": family + ":italic", "mathtext.bf": family + ":bold",
                         "mathtext.cal": family,
                         "axes.labelcolor": INK, "svg.fonttype": "path",
                         "svg.hashsalt": "unwrap-locust-density-v1", "figure.facecolor": "white"}):
        fig = plt.figure(figsize=(8.8, 7.1))
        photos = []
        heading(fig, "a", "Real locust", .025, .973, .135)
        heading(fig, "b", "Virtual locust", .270, .973, .365)
        heading(fig, "c", "Immersive arena", .495, .973, .74)
        photos.append(photo(fig, "real_locust", [.025, .725, .22, .20], crop=[.20, 0, .80, 1]))
        # Original front/lateral/rear views, with generous space around limbs.
        for name, y in (("virtual_locust_front", .858), ("virtual_locust_lateral", .791), ("virtual_locust_rear", .725)):
            photos.append(photo(fig, name, [.270, y, .195, .067],
                                crop=[.18, .13, .82, .89], fit="contain"))
        photos.append(photo(fig, "vr_arena", [.495, .725, .48, .20]))
        heading(fig, "d", "Swarm stimuli", .025, .676, .135)
        heading(fig, "e", "Density & order", .270, .676, .372)
        heading(fig, "f", "Model comparison", .490, .676, .623)
        heading(fig, "g", "Fixed effects", .765, .676, .88)
        # Match the -90-degree image orientation in the source PDF, so the
        # displayed virtual agents and the marching arrow point together.
        photos.append(photo(fig, "real_swarm", [.025, .390, .22, .215], rotation=-90))
        photos.append(photo(fig, "virtual_swarm", [.025, .125, .22, .215], rotation=-90))
        fig.text(.135, .622, "Real swarm", ha="center", fontsize=7.5)
        fig.text(.135, .356, "Virtual swarm", ha="center", fontsize=7.5)
        arrow = fig.add_axes([.075, .071, .12, .026]); arrow.set_axis_off()
        arrow.annotate("", xy=(1, .6), xytext=(0, .6), arrowprops={"arrowstyle": "-|>", "lw": 1.2, "color": INK})
        fig.text(.135, .046, "Marching direction", ha="center", fontsize=7.5)
        stimulus_schematic(fig)
        # Forest plots: intervals first and markers second; no gridlines.
        ax = fig.add_axes([.574, .435, .135, .17]); style_axis(ax)
        rows = selection["models"]
        for y, row in enumerate(rows[::-1]):
            color = "#1f77b4" if row["has_order"] else "#ff7f0e"
            ax.plot([row["delta"] - row["se"], row["delta"] + row["se"]], [y, y], color=color, lw=1.1)
            ax.scatter(row["delta"], y, s=22, color=color, zorder=3,
                       marker="o" if row["has_order"] else "s")
        ax.axvline(0, lw=.7, color=MUTED, ls=(0, (3, 3)))
        ax.set_yticks(range(5), [r["name"].replace("Density", "D").replace("Order", "O") for r in rows[::-1]])
        ax.set_ylim(-.5, 4.6); ax.set_xlim(-96, 6); ax.set_xticks([-80, -40, 0])
        ax.set_xlabel("ELPD difference", labelpad=5)
        ax.legend(handles=[Line2D([], [], color=c, marker=m, ls="", label=t, ms=4)
                           for c, m, t in (("#1f77b4", "o", "Order included"), ("#ff7f0e", "s", "Order excluded"))],
                  loc="lower left", bbox_to_anchor=(-.04, 1.03), frameon=False,
                  borderpad=0, fontsize=6.5, handletextpad=.4)
        ax = fig.add_axes([.835, .435, .140, .17]); style_axis(ax)
        for term, row in enumerate(stats["coefficients"]):
            for j, (color, marker) in enumerate(((MU_COLOR, "D"), (KAPPA_COLOR, "^"))):
                y = 3 - term + (.10 if j == 0 else -.10)
                mean, low, high = row[j]
                ax.plot([low, high], [y, y], color=color, lw=1.1)
                ax.scatter(mean, y, color=color, marker=marker, s=22, zorder=3)
        ax.axvline(0, lw=.7, color=MUTED, ls=(0, (3, 3)))
        ax.set_yticks(range(4), ["D × O", "Order", "Density", "Intercept"])
        ax.set_ylim(-.5, 3.6); ax.set_xlabel("Coefficient value", labelpad=5)
        ax.legend(handles=[Line2D([], [], color=c, marker=m, ls="", label=t, ms=4)
                           for c, m, t in ((MU_COLOR, "D", "Direction μ"), (KAPPA_COLOR, "^", "Concentration κ"))],
                  loc="lower left", bbox_to_anchor=(-.05, 1.03), frameon=False,
                  borderpad=0, fontsize=6.5, handletextpad=.4)
        heading(fig, "h", "Response to swarm motion", .490, .344, .75)
        for key, grid, bounds, levels, cmap, xlabel in (
            ("over_order", "order_grid", [.535, .105, .17, .155], DENSITY_LEVELS, plt.cm.Reds, "Order"),
            ("over_density", "density_grid", [.805, .105, .17, .155], ORDER_LEVELS, plt.cm.Purples, "Density (m$^{-2}$)"),
        ):
            ax = fig.add_axes(bounds); style_axis(ax)
            for i, (level, curve) in enumerate(zip(levels, stats[key])):
                color = cmap(.45 + .10 * i)
                ax.fill_between(stats[grid], curve["low"], curve["high"], color=color, alpha=.15, lw=0)
                ax.plot(stats[grid], curve["mean"], color=color, lw=1.25)
            ax.axhline(1 / (2 * np.pi), color=MUTED, lw=.7, ls=(0, (3, 3)))
            ax.set_ylim(0, 1); ax.set_yticks([0, .25, .5, .75, 1], ["0", ".25", ".5", ".75", "1"])
            ax.set_xlabel(xlabel, labelpad=5)
            if key == "over_order":
                ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1]); ax.set_ylabel("Density at 0° (rad$^{-1}$)", labelpad=5)
            else:
                ax.set_xscale("log", base=2); ax.set_xlim(1, 64)
                ax.set_xticks([1, 4, 16, 64], ["1", "4", "16", "64"])
            # Compact colour key above each plot; all six levels are explicit.
            handles = [Line2D([], [], color=cmap(.45 + .10*i), lw=2, label=f"{v:g}") for i,v in enumerate(levels)]
            ax.legend(handles=handles, title="Density (m$^{-2}$)" if key=="over_order" else "Order",
                      loc="lower left", bbox_to_anchor=(-.05, 1.01), ncol=3, frameon=False,
                      fontsize=6.3, title_fontsize=6.8, handlelength=.8,
                      columnspacing=.7, handletextpad=.35, labelspacing=.25, borderpad=0)
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg", metadata={"Date": None, "Creator": "unwrap locust_density_figure.py"})
        plt.close(fig)
    return embed_photos(buffer.getvalue(), photos)


def embed_photos(svg, placements):
    """Use original JPEG streams inside vector SVG frames and clipping paths."""
    ET.register_namespace("", SVG)
    root = ET.fromstring(svg)
    title = ET.Element(f"{{{SVG}}}title", id="figure-title"); title.text = "Locust density and order"
    desc = ET.Element(f"{{{SVG}}}desc", id="figure-description"); desc.text = CAPTION
    root.insert(0, title); root.insert(1, desc)
    root.set("role", "img"); root.set("aria-labelledby", "figure-title figure-description")
    manifest = json.loads((ASSETS / "images.json").read_text())
    width, height = 8.8 * 72, 7.1 * 72
    for i, placement in enumerate(placements):
        name = placement["name"] + ".jpeg"
        path = ASSETS / name
        meta = manifest["images"][name]
        if sha256(path) != meta["sha256"]:
            raise ValueError(f"Altered source image: {path}")
        x, y, w, h = placement["bounds"]
        x, y, w, h = x * width, (1 - y - h) * height, w * width, h * height
        crop = placement["crop"] or [0, 0, 1, 1]
        rotation = placement["rotation"]
        source_w, source_h = meta["width"], meta["height"]
        if rotation == -90:
            source_w, source_h = source_h, source_w
        # Crop to an explicit source rectangle, then fill the destination;
        # preserve the source aspect ratio and clip any extra overflow.
        l, t, r, b = crop
        fit = min if placement["fit"] == "contain" else max
        scale = fit(w / (source_w * (r-l)), h / (source_h * (b-t)))
        iw, ih = source_w * scale, source_h * scale
        ix, iy = x + w/2 - (l+r)/2*iw, y + h/2 - (t+b)/2*ih
        clip = ET.SubElement(root, f"{{{SVG}}}clipPath", id=f"photo-clip-{i}")
        ET.SubElement(clip, f"{{{SVG}}}rect", x=str(x), y=str(y), width=str(w), height=str(h))
        group = ET.SubElement(root, f"{{{SVG}}}g", {"clip-path": f"url(#photo-clip-{i})"})
        ET.SubElement(group, f"{{{SVG}}}rect", x=str(x), y=str(y), width=str(w), height=str(h), fill="#e5e5e5")
        attrs = {"x": str(ix), "y": str(iy), "width": str(iw), "height": str(ih)}
        if rotation == -90:
            attrs = {"x": "0", "y": "0", "width": str(ih), "height": str(iw),
                     "transform": f"translate({ix} {iy + ih}) rotate(-90)"}
        ET.SubElement(group, f"{{{SVG}}}image", {**attrs,
            "href": "data:image/jpeg;base64," + base64.b64encode(path.read_bytes()).decode(),
            "preserveAspectRatio": "xMidYMid meet"})
    return "\n".join(line.rstrip() for line in ET.tostring(root, encoding="unicode").splitlines()) + "\n"


def write_outputs(output, stats, selection, metadata, *, family="Arial"):
    import cairosvg
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    stem = output / "locust_density_publication"
    svg = draw_figure(stats, selection, family=family)
    stem.with_suffix(".svg").write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(stem.with_suffix(".pdf")))
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(stem.with_suffix(".png")), output_width=2200)
    stem.with_suffix(".html").write_text(f'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Locust density and order</title>
<style>body{{margin:0;background:#f5f5f4;color:#262626;font:15px/1.65 Arial,sans-serif}}
main{{max-width:1180px;margin:40px auto;padding:0 24px}}figure{{margin:0}}.figure{{background:white;padding:20px;border:1px solid #e6e6e6}}svg{{display:block;width:100%;height:auto}}figcaption{{max-width:96ch;margin:24px auto}}a{{color:#665017}}</style>
<main><nav><a href="locust_density_publication.svg">SVG</a> · <a href="locust_density_publication.pdf">PDF</a></nav>
<figure><div class="figure">{svg}</div><figcaption>{CAPTION}</figcaption></figure></main></html>''')
    metadata["summaries"] = {"coefficients": stats["coefficients"].tolist(),
        "density_at_zero_order_endpoints": [dict(density=float(d), geometric_mean=curve["mean"][[0,-1]].tolist(),
            low=curve["low"][[0,-1]].tolist(), high=curve["high"][[0,-1]].tolist())
            for d,curve in zip(DENSITY_LEVELS,stats["over_order"])]}
    metadata["outputs_sha256"] = {s: sha256(stem.with_suffix(s)) for s in (".svg", ".pdf", ".png", ".html")}
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2) + "\n")
    return stem


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=ROOT / ".cache/locust-figure" / DATA_REVISION)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/locust-density")
    parser.add_argument("--font-dir", type=Path)
    parser.add_argument("--font-family", default="Arial")
    args = parser.parse_args(argv)
    fonts = resolve_fonts(args.font_dir, args.font_family)
    root = prepare_inputs(args.data_root, download=args.download)
    path = root / "fit_cache/locust_density_hier.nc"
    with xr.open_dataset(path, group="posterior", engine="h5netcdf") as ds:
        post = ds.load()
    with xr.open_dataset(path, group="sample_stats", engine="h5netcdf") as ds:
        divergences = int(ds.diverging.sum())
    if divergences:
        raise ValueError(f"Cached fit has {divergences} divergences; review before plotting")
    selection = json.loads((ASSETS / "selection.json").read_text())
    sources = [Path(__file__), ROOT / "tools/prepare_locust_figure_assets.py",
               ROOT / "unwrap/link/numpy.py", ROOT / "unwrap/circular.py",
               ROOT / "unwrap/bessel.py", *sorted(ASSETS.iterdir())]
    metadata = {"data_repository": "jgraving/unwrap-data", "data_revision": DATA_REVISION,
        "inputs_sha256": INPUT_HASHES, "code_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "source_sha256": {p.relative_to(ROOT).as_posix(): sha256(p) for p in sources if p.is_file()},
        "python": platform.python_version(), "font_family": args.font_family, "fonts": fonts,
        "packages": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "matplotlib", "xarray", "h5netcdf", "h5py", "cairosvg")},
        "posterior_chains": post.sizes["chain"], "posterior_draws_per_chain": post.sizes["draw"],
        "divergences": divergences, "selection_source": selection,
        "prediction_basis": "Original scaled coordinates of stored Bambi coefficients; concentration offset zero",
        "prediction_summary": "exp(mean(log density at zero)); exponentiated 90% log-density HDI",
        "coefficient_hdi_probability": .95, "schematic_seed": 123, "caption": CAPTION}
    stem = write_outputs(args.output_dir, statistics(post), selection, metadata, family=args.font_family)
    print(f"Wrote SVG, PDF, PNG, HTML and provenance JSON: {stem}")


if __name__ == "__main__":
    main()
