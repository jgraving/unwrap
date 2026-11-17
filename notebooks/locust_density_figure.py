"""Rebuild the locust density/order publication figure without sampling.

Run from the repository root; see LOCUST_DENSITY_FIGURE.md for provenance,
the corrected predictor coordinates, and source-image extraction.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import io
import json
from pathlib import Path
import platform
import re
import subprocess
import urllib.request
import xml.etree.ElementTree as ET

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
import numpy as np
from scipy.special import i0e
import xarray as xr

from unwrap.circular import circular_mean
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
    for side, spine in ax.spines.items():
        spine.set_visible(side in ("left", "bottom"))
        spine.set_linewidth(.65)
        spine.set_color(INK)
    ax.tick_params(length=3, width=.65, pad=3, colors=INK, labelsize=9)
    ax.grid(False)


def draw_figure(stats, selection, *, family="Arial"):
    """Polish the manuscript's 750 × 525-point multipanel layout.

    Source photo transforms and the calibrated scale bar are retained.
    The presentation pass lightens frames, updates typography, and redraws
    f–h from the same numerical inputs. Units remain single math expressions.
    """
    width, height = 750, 525
    def axes_at(fig, x, y, w, h):
        return fig.add_axes([x/width, (height-y-h)/height, w/width, h/height])

    with mpl.rc_context({"font.family": family, "font.size": 9,
                         "text.color": INK, "axes.labelcolor": INK,
                         "mathtext.fontset": "custom", "mathtext.rm": family,
                         "mathtext.it": family + ":italic", "mathtext.bf": family + ":bold",
                         "mathtext.cal": family, "svg.fonttype": "path",
                         "svg.hashsalt": "unwrap-locust-density-polish-v1"}):
        fig = plt.figure(figsize=(width/72, height/72))
        fig.patch.set_visible(False)
        # Retain the original title wording and positions with a common weight.
        for text, x, y in (("Model selection", 471.45, 198.81),
                           ("Population-level\ncoefficients", 663.04, 193.39),
                           ("Population-level summary", 547.97, 377.55)):
            fig.text(x/width, 1-y/height, text, ha="center", va="top",
                     fontsize=11, weight="bold", linespacing=1)
        fig.text(245.4/width, 1-225/height, r"$(\mathrm{m}^{-2})$",
                 ha="center", va="baseline", fontsize=11)
        ax = axes_at(fig, 420.534, 220.492, 107.159, 116.196)
        style_axis(ax)
        rows = selection["models"]
        for y, row in enumerate(rows[::-1]):
            color = "#1f77b4" if row["has_order"] else "#ff7f0e"
            ax.errorbar(row["delta"], y, xerr=row["se"], color=color,
                        fmt="o" if row["has_order"] else "s", ms=5.7,
                        capsize=3, elinewidth=1.1, capthick=1)
        ax.axvline(0, lw=.8, color=MUTED, ls=(0, (4, 3)))
        ax.set_yticks(range(5), [r["name"] for r in rows[::-1]])
        ax.set_ylim(-.2, 4.2); ax.set_xlim(-95, 5); ax.set_xticks([-50, 0])
        ax.set_xlabel("ELPD difference", labelpad=6)
        ax.legend(handles=[Line2D([], [], color=c, marker=m, label=t, ms=5.5)
                           for c, m, t in (("#1f77b4", "o", "Order included"),
                                           ("#ff7f0e", "s", "Order excluded"))],
                  loc="upper left", frameon=False, fontsize=9, borderpad=.25,
                  handlelength=.8, handletextpad=.35, labelspacing=.15)
        ax = axes_at(fig, 607.833, 220.492, 107.159, 116.196)
        style_axis(ax)
        for term, row in enumerate(stats["coefficients"]):
            for j, (color, marker) in enumerate(((MU_COLOR, "D"), (KAPPA_COLOR, "^"))):
                y = 3-term + (.075 if j == 0 else -.075)
                mean, low, high = row[j]
                midpoint = (low+high)/2
                ax.errorbar(midpoint, y, xerr=(high-low)/2, color=color,
                            fmt="none", capsize=3, elinewidth=1.1, capthick=1)
                ax.plot(mean, y, marker, ms=6, color=color)
        ax.axvline(0, lw=.8, color=MUTED, ls=(0, (4, 3)))
        ax.set_yticks(range(4), ["Density × Order", "Order", "Density", "Intercept"])
        ax.set_ylim(-.25, 3.25); ax.set_xlim(-.35, 3.02); ax.set_xticks([0, 2])
        ax.set_xlabel("Coefficient value", labelpad=6)
        ax.legend(handles=[Line2D([], [], color=c, marker=m, label=t, ms=6)
                           for c, m, t in ((MU_COLOR, "D", "μ"), (KAPPA_COLOR, "^", "κ"))],
                  loc="upper right", frameon=False, fontsize=9, borderpad=.25,
                  handlelength=.8, handletextpad=.35, labelspacing=.15)
        for key, grid, x, levels, cmap, xlabel in (
            ("over_order", "order_grid", 405.066, DENSITY_LEVELS, plt.cm.Reds, "Order"),
            ("over_density", "density_grid", 587.715, ORDER_LEVELS, plt.cm.Purples, "Density"),
        ):
            ax = axes_at(fig, x, 401.301, 110.682, 89.75)
            style_axis(ax)
            handles = []
            for i in range(len(levels)-1, -1, -1):
                curve = stats[key][i]
                color = cmap(.45 + .10*i)
                ax.fill_between(stats[grid], curve["low"], curve["high"], color=color, alpha=.16, lw=0)
                line, = ax.plot(stats[grid], curve["mean"], color=color, lw=1.5)
                handles.append(line)
            ax.axhline(1/(2*np.pi), color=MUTED, lw=.8, ls=(0, (4, 3)))
            ax.set_ylim(0, 1); ax.set_yticks(np.linspace(0, 1, 6))
            ax.set_xlabel(xlabel, labelpad=6)
            if key == "over_order":
                ax.set_xlim(0, 1); ax.set_xticks(np.linspace(0, 1, 6))
                ax.set_xticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"])
                ax.set_ylabel(r"Predictive density $(\mathrm{rad}^{-1})$", labelpad=4)
            else:
                ax.set_xscale("log", base=2); ax.set_xlim(1, 64)
                ax.set_xticks([1, 2, 4, 8, 16, 32, 64], ["1", "2", "4", "8", "16", "32", "64"])
            labels = [f"{v:g}" if key == "over_order" else f"{v:.1f}" for v in levels[::-1]]
            ax.legend(handles=handles, labels=labels,
                      title="Density" if key=="over_order" else "Order",
                      loc="center left", bbox_to_anchor=(1.035, .5), ncol=1,
                      frameon=False, fontsize=9, title_fontsize=9,
                      handlelength=1.2, handletextpad=.4, labelspacing=.15, borderpad=.3)
        buffer = io.StringIO()
        fig.savefig(buffer, format="svg", metadata={"Date": None, "Creator": "unwrap locust_density_figure.py"})
        plt.close(fig)
    return embed_original_panels(buffer.getvalue())


def polish_original_panels(root):
    """Style selected vector elements without moving photos or the scale bar.

    Operate on a parsed copy; the checksum-verified reproduction asset stays
    intact. The source SVG retains its original PDF coordinates and clips.
    """
    parents = {child: parent for parent in root.iter() for child in parent}
    anchors = [(28.97936, 25.457), (237.8304, 25.739808), (377.0035, 24.89151),
               (28.45983, 209.1947), (189.6003, 209.1947), (351.6703, 209.1947),
               (557.7246, 209.1947), (349.0309, 392.63349)]
    changed_letters = removed_titles = schematic_borders = photo_borders = 0
    for node in list(root.iter()):
        transform = node.get("transform", "")
        match = re.fullmatch(r"matrix\(([^)]+)\)", transform)
        values = [float(v) for v in re.split(r"[ ,]+", match[1])] if match else []
        if node.tag == f"{{{SVG}}}use" and len(values) == 6:
            size, _, _, _, x, y = values
            if abs(size - 19) < .001:
                ax, ay = min(anchors, key=lambda a: abs(a[0]-x) + 10*abs(a[1]-y))
                scale = 16/19
                node.set("transform", f"matrix(16,0,0,-16,{ax+(x-ax)*scale},{ay-1.5})")
                changed_letters += 1
            elif abs(size - 11) < .001 and x > 420 and (198 < y < 215 or 385 < y < 390):
                parents[node].remove(node)
                removed_titles += 1
        if node.tag != f"{{{SVG}}}path" or node.get("stroke") != "#000000":
            continue
        path = node.get("d", "")
        if transform.startswith("matrix(1.651933,") and node.get("stroke-width") == ".8" \
                and re.fullmatch(r"M[\d.]+ [\d.]+[VH][\d.]+", path):
            node.set("stroke-width", ".4")
            node.set("stroke", "#858585")
            schematic_borders += 1
        elif node.get("stroke-width") == "1.5" and path == "M0 0V119.125H119.4375V0H0ZM0 0":
            node.set("stroke-width", ".75")
            node.set("stroke", "#858585")
            photo_borders += 1
    if (changed_letters, schematic_borders, photo_borders) != (16, 48, 2) or removed_titles < 50:
        raise ValueError("Original SVG structure changed; review the targeted presentation edits")


def embed_original_panels(svg):
    """Compose preserved a–e artwork underneath newly computed f–h panels."""
    ET.register_namespace("", SVG)
    root = ET.fromstring(svg)
    manifest = json.loads((ASSETS / "images.json").read_text())
    info = manifest["original_panels"]
    path = ASSETS / info["file"]
    if sha256(path) != info["sha256"]:
        raise ValueError("Original panel artwork checksum mismatch")
    original = ET.parse(path).getroot()
    polish_original_panels(original)
    original.set("x", "0"); original.set("y", "0")
    title = ET.Element(f"{{{SVG}}}title", id="figure-title"); title.text = "Locust density and order"
    desc = ET.Element(f"{{{SVG}}}desc", id="figure-description"); desc.text = CAPTION
    background = ET.Element(f"{{{SVG}}}rect", x="0", y="0", width="750", height="525", fill="white")
    for index, item in enumerate((title, desc, background, original)):
        root.insert(index, item)
    root.set("role", "img"); root.set("aria-labelledby", "figure-title figure-description")
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
        "coefficient_hdi_probability": .95, "layout": "Polished original 750 × 525-point manuscript layout",
        "illustration_basis": "Original PDF artwork for panels a–e, including calibrated scale bar",
        "presentation": {"panel_letter_points": 16, "title_points": 11, "title_weight": "bold",
                         "statistical_label_points": 9, "prediction_line_points": 1.5,
                         "prediction_band_alpha": .16, "legend_columns": 1, "legend_frames": False},
        "caption": CAPTION}
    stem = write_outputs(args.output_dir, statistics(post), selection, metadata, family=args.font_family)
    print(f"Wrote SVG, PDF, PNG, HTML and provenance JSON: {stem}")


if __name__ == "__main__":
    main()
