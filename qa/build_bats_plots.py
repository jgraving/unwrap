"""Build goldshtein_bats_plots.py by slicing the annotated cells verbatim.

Every figure becomes a module function that takes the fit results and the
modeling table, and saves through the house `_save`/`unwrap.plots` machinery
(100 dpi inline in the notebook, 300 dpi file via an explicit save_path).
Cell-local helpers are deduped: identical copies are dropped (the module
imports the shared one from goldshtein_bats_posthoc), divergent copies are
renamed with a suffix naming the figure family they belong to.
"""
import re

DUMP = "/tmp/opencode/bats_cells.txt"
txt = open(DUMP).read()

cells = {}
for m in re.finditer(
    r"########## (c\d+) (code|markdown) ##########\n(.*?)(?=\n########## |\Z)", txt, re.S
):
    if m.group(2) == "code":
        cells[f"c{int(m.group(1)[1:]):02d}"] = m.group(3)


def L(i):
    return cells[f"c{i:02d}"].split("\n")


def sub1(blob, old, new, what):
    c = blob.count(old)
    if c == 0 and old.endswith("\n"):
        # slice boundary: def ends at the cell text, no trailing newline
        old, new = old.rstrip("\n"), new.rstrip("\n")
        c = blob.count(old)
    if c != 1:
        import difflib

        first = old.split("\n")[0]
        i = blob.find(first)
        region = blob[i : i + len(old) + 120] if i >= 0 else "(first line absent)"
        d = "\n".join(
            difflib.unified_diff(old.split("\n"), region.split("\n"), "pinned", "actual", lineterm="")
        )
        raise AssertionError(f"{what}: found {c} (expected 1); near-diff:\n{d}")
    return blob.replace(old, new)


def subn(blob, old, new, n, what):
    c = blob.count(old)
    assert c == n, f"{what}: found {c} (expected {n}): {old[:120]!r}"
    return blob.replace(old, new)


def drop_line(blob, marker, what):
    ls = blob.split("\n")
    hits = [k for k, l in enumerate(ls) if l.strip() and marker in l]
    assert len(hits) == 1, f"{what}: {len(hits)} lines contain {marker!r}"
    del ls[hits[0]]
    return "\n".join(ls)


def topdef(i, name):
    """(start, signature-close-inclusive) of a top-level def."""
    ls = L(i)
    s = next((k for k, l in enumerate(ls) if l.startswith(f"def {name}(")), None)
    assert s is not None, f"c{i:02d}: def {name} not found"
    if ls[s].rstrip().endswith("):"):
        return s, s  # single-line signature
    bal = 0
    e = s
    for j in range(s, len(ls)):
        bal += ls[j].count("(") - ls[j].count(")")
        if bal <= 0 and j > s:
            e = j
            break
    return s, e


def def_range(i, name):
    """(start, end-exclusive) of a top-level def: body runs to the next
    top-level line (trailing blanks/comments belong to the following block)."""
    ls = L(i)
    s, e = topdef(i, name)
    be = len(ls)
    for j in range(e + 1, len(ls)):
        l = ls[j]
        if l.strip() and not l.startswith((" ", "\t")):
            be = j
            break
    while be > s + 1 and (not ls[be - 1].strip() or ls[be - 1].lstrip().startswith("#")):
        be -= 1
    return s, be


def body(i, a, b):
    """Inclusive 1-based line slice of a cell."""
    ls = L(i)
    src = "\n".join(ls[a - 1 : b])
    assert src.strip(), f"c{i:02d} {a}..{b} empty"
    return src


def indent(src, pad="    "):
    return "\n".join(pad + l if l.strip() else "" for l in src.split("\n"))


def norm(s):
    """Canonical form for equality tests (docstrings + blank lines stripped)."""
    s = re.sub(r'"""[\s\S]*?"""', "", s)
    return "\n".join(l.rstrip() for l in s.split("\n") if l.strip()).strip()


OUT = []
add = OUT.append

# ---------------------------------------------------------------- header
hdr = open("/tmp/opencode/plots_header.txt").read()
hdr = sub1(
    hdr,
    "from matplotlib.lines import Line2D\n",
    "from matplotlib.collections import LineCollection\nfrom matplotlib.lines import Line2D\n",
    "header LineCollection",
)
hdr = sub1(
    hdr,
    "    log_i0_stable,\n    log_vonmises_pdf,\n)",
    "    log_i0_stable,\n    log_vonmises_pdf,\n    build_res_treat_from_vonmises_idata,\n    summarize_band,\n    format_polar,\n    format_treatment_title,\n)",
    "header build_res import",
)
hdr += (
    "\n"
    "N_COLS = 8  # grid columns of the individual PPC polar grid\n"
    "\n"
    "#: radial constants of the 2x2 mixture-polar family "
    "(plot_mixture_polar_pdfs and its _polar_format_mix helper)\n"
    "EDGE_RADIUS = 1.0\n"
    "AMP = 1.0\n"
    "OUTER_LW = 0.8\n"
)
add(hdr.rstrip("\n"))

# ============================================================ map toolkit (c25)
toolkit = [
    "_nice_scale_length_m",
    "_pick_order_cols",
    "_plot_gdf_capture_artists",
    "_set_rasterized",
    "_style_trajectory_artists",
    "add_scale_bar_dual",
    "add_compass",
    "_add_scale_and_compass",
    "_ensure_dest_cols",
    "_dest_label_from_flags",
    "_make_segment_and_gap_lines_from_df",
]
toolkit = ["hide_axes"] + toolkit
for name in toolkit:
    s, e = def_range(25, name)
    add("\n".join(L(25)[s:e]))

# ============================================================ c27 -> plot_destination_panels
s, e = def_range(27, "plot_rediscretized_destination_panels")
fn27 = "\n".join(L(27)[s:e])
fn27 = sub1(fn27, "def plot_rediscretized_destination_panels(", "def plot_destination_panels(", "c27 rename")
fn27 = sub1(
    fn27,
    '    out_pdf="rediscretized_destination_panels.pdf",\n    out_svg="rediscretized_destination_panels.svg",',
    "",
    "c27 drop out paths",
)
fn27 = sub1(
    fn27,
    "    save_dpi: int = 300,                     # affects raster layers in PDF/SVG export\n",
    "",
    "c27 drop save_dpi",
)
fn27 = sub1(
    fn27,
    "    show_treatment_legend: bool = True,\n):",
    "    show_treatment_legend: bool = True,\n    save: bool = True,\n    save_path=None,\n):",
    "c27 save kwargs",
)
fn27 = sub1(
    fn27,
    """        os.makedirs(os.path.dirname(out_pdf), exist_ok=True)
        os.makedirs(os.path.dirname(out_svg), exist_ok=True)

        # Note: for raster tiles, real detail mainly comes from basemap_zoom / retina.
        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02, dpi=save_dpi)
        fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.02, dpi=save_dpi)

        plt.show()
        plt.close(fig)""",
    "        _save(fig, save, save_path)\n        return fig",
    "c27 tail",
)
add(fn27)

# ============================================================ c25 -> plot_destination_panels_compact
s, e = def_range(25, "plot_rediscretized_destination_panels")
fn25 = "\n".join(L(25)[s:e])
fn25 = sub1(fn25, "def plot_rediscretized_destination_panels(", "def plot_destination_panels_compact(", "c25 rename")
fn25 = sub1(
    fn25,
    '    out_pdf="../figures/bats/rediscretized_destination_panels_compact.pdf",\n    out_svg="../figures/bats/rediscretized_destination_panels_compact.svg",\n',
    "",
    "c25 drop out paths",
)
fn25 = sub1(
    fn25,
    """    # selective rasterization controls
    rasterize_trajectories_pdf=True,
    rasterize_trajectories_svg=False,   # keep SVG trajectories vector by default (nicest)
""",
    "",
    "c25 drop rasterize flags",
)
fn25 = sub1(
    fn25,
    "    svg_raster_dpi=600,                 # only matters if rasterize_trajectories_svg=True\n):",
    "    svg_raster_dpi=600,  # raster dpi of the single saved PDF export\n    save: bool = True,\n    save_path=None,\n):",
    "c25 save kwargs",
)
fn25 = sub1(
    fn25,
    """    # --------------------------------------------------------
    # SAVE: selectively rasterize ONLY trajectories per export
    # --------------------------------------------------------
    if out_pdf:
        _set_rasterized(trajectory_artists, rasterize_trajectories_pdf)
        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.01, dpi=pdf_raster_dpi)

    if out_svg:
        _set_rasterized(trajectory_artists, rasterize_trajectories_svg)
        fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.01, dpi=svg_raster_dpi)

    plt.show()
    return fig, axes""",
    "    _rasterize = str(save_path).endswith(\".pdf\")\n"
    "    _set_rasterized(trajectory_artists, _rasterize)\n"
    "    _save(fig, save, save_path, dpi=pdf_raster_dpi if _rasterize else svg_raster_dpi)\n"
    "    return fig, axes",
    "c25 tail",
)
add(fn25)

# ============================================================ c29 flow fields
# c29 duplicates of the c25 toolkit must either match (drop) or be renamed
for name in ("_ensure_dest_cols", "_dest_label_from_flags", "_pick_order_cols"):
    a = "\n".join(L(29)[def_range(29, name)[0] : def_range(29, name)[1]])
    b = "\n".join(L(25)[def_range(25, name)[0] : def_range(25, name)[1]])
    assert norm(a) == norm(b), f"c29 duplicate {name} differs from c25"
s, e = def_range(29, "_nice_scale_length_m")
nice29 = "\n".join(L(29)[s:e])
nice29 = sub1(nice29, "def _nice_scale_length_m(", "def _nice_scale_length_m_flow(", "c29 nice rename")
add(nice29)
for name in (
    "_format_scale_label",
    "_fallback_add_scale_and_compass",
    "_safe_normalize",
    "_make_lines_from_df",
    "lake_outer_polygon_largest",
    "polygon_mask_on_grid",
    "_adaptive_grid_res",
    "solve_laplace_dirichlet",
    "harmonic_field_on_extent",
    "smoothed_distance_field_on_extent",
    "roost_field_on_extent",
):
    s, e = def_range(29, name)
    blk = "\n".join(L(29)[s:e])
    blk = subn(blk, "_nice_scale_length_m(", "_nice_scale_length_m_flow(", blk.count("_nice_scale_length_m("), f"c29 {name} nice call")
    add(blk)

s, e = def_range(29, "plot_destination_flow_fields_panels")
fn29 = "\n".join(L(29)[s:e])
fn29 = sub1(fn29, "def plot_destination_flow_fields_panels(", "def plot_destination_flow_fields(", "c29 rename")
fn29 = drop_line(fn29, 'out_pdf="destination_flow_fields.pdf"', "c29 out_pdf")
fn29 = drop_line(fn29, 'out_svg="destination_flow_fields.svg"', "c29 out_svg")
fn29 = sub1(
    fn29,
    "    quiver_minshaft: float = 2.2,\n):",
    "    quiver_minshaft: float = 2.2,\n    save: bool = True,\n    save_path=None,\n):",
    "c29 save kwargs",
)
fn29 = sub1(
    fn29,
    """    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.02, dpi=300)
    fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.show()
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c29 tail",
)
fn29 = subn(fn29, "_nice_scale_length_m(", "_nice_scale_length_m_flow(", fn29.count("_nice_scale_length_m("), "c29 flow nice call")
add(fn29)

# ============================================================ c31 per-tag polar KDE
s, e = def_range(31, "plot_polar_vonmises_kde_all")
fn31 = "\n".join(L(31)[s:e])
fn31 = sub1(fn31, "def plot_polar_vonmises_kde_all(", "def plot_per_tag_polar_kde(", "c31 rename")
fn31 = sub1(
    fn31,
    "    title_fs=9,\n):",
    "    title_fs=9,\n    save: bool = True,\n    save_path=None,\n):",
    "c31 save kwargs",
)
fn31 = sub1(
    fn31,
    "    plt.tight_layout()\n    plt.show()\n",
    "    plt.tight_layout()\n    _save(fig, save, save_path)\n    return fig\n",
    "c31 tail",
)
fn31 = sub1(fn31, 'cmap = plt.get_cmap("tab10")', 'cmap = mpl.colormaps.get_cmap("tab10")', "c31 cmap")
add(fn31)

# ============================================================ c50 coefficient plot
c50_body = body(50, 27, 244)
c50_body = sub1(c50_body, "h = az.hdi(centered, hdi_prob=hdi_prob, skipna=True)", "h = az.hdi(centered, prob=hdi_prob, skipna=True)", "c50 az.hdi a")
c50_body = sub1(c50_body, "h = az.hdi(x, hdi_prob=hdi_prob, skipna=True)", "h = az.hdi(x, prob=hdi_prob, skipna=True)", "c50 az.hdi b")
c50_body = sub1(c50_body, 'cmap = mpl.cm.get_cmap("tab10")', 'cmap = mpl.colormaps.get_cmap("tab10")', "c50 cmap")
c50_body = sub1(
    c50_body,
    "plt.figure(figsize=(1.8, 0.40 * len(terms_order)), dpi=300)",
    "fig = plt.figure(figsize=(1.8, 0.40 * len(terms_order)), dpi=300)",
    "c50 fig",
)
c50_body = sub1(
    c50_body,
    'plt.savefig("../figures/bats/treatment_coefficients_mu_kappa.pdf", dpi=300, bbox_inches="tight", pad_inches=0.01)\nplt.show()\n',
    "_save(fig, save, save_path)\nreturn fig\n",
    "c50 save",
)
uses_df = "df_track" in c50_body
add(
    "def plot_treatment_coefficients_mu_kappa(idata_treatment"
    + (", df_track" if uses_df else "")
    + ", *, save=True, save_path=None):\n"
    "    \"\"\"Population mu/kappa coefficients of the Bambi von Mises model "
    "(intercept + treatment effects) with circular/linear mean + 95% HDI.\"\"\"\n"
    + indent(c50_body)
)

# ============================================================ c52 predictive density at 0
c52_body = body(52, 40, 266)
c52_body = sub1(
    c52_body,
    "hdi  = np.array([az.hdi(d, hdi_prob=HDI_PROB) for d in dens_draws_by_row], float)  # (K, 2)",
    "hdi  = np.array([az.hdi(d, prob=HDI_PROB) for d in dens_draws_by_row], float)  # (K, 2)",
    "c52 az.hdi",
)
c52_body = sub1(
    c52_body,
    """plt.tight_layout()
plt.savefig(
    "../figures/bats/treatment_pred_density_at0_mu_kappa.pdf",
    dpi=DPI,
    bbox_inches="tight",
    pad_inches=0.01
)
plt.show()
""",
    """plt.tight_layout()
_save(fig, save, save_path)
return fig
""",
    "c52 save",
)
add(
    "def plot_treatment_pred_density_at0_mu_kappa(idata_treatment, df_track, *, save=True, save_path=None):\n"
    "    \"\"\"Predictive density at 0 deg of the Bambi von Mises model for the population "
    "intercept and each treatment shift (mean + 95% HDI, uniform reference line).\"\"\"\n"
    + indent(c52_body)
)

# ============================================================ PSIS-LOO set
# c59 -> plot_total_elpd
s, e = def_range(59, "plot_total_elpd_from_summary")
fn59 = "\n".join(L(59)[s:e])
fn59 = sub1(fn59, "def plot_total_elpd_from_summary(", "def plot_total_elpd(", "c59 rename")
fn59 = sub1(fn59, "    *,\n    savepath=None,", "    *,\n    save=True,\n    save_path=None,", "c59 sig")
fn59 = sub1(
    fn59,
    """    if savepath:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight", pad_inches=0.01)
    plt.show()
    return fig, ax
""",
    """    _save(fig, save, save_path)
    return fig, ax
""",
    "c59 tail",
)
add(fn59)

# c61 -> plot_delta_total_and_by_treatment
s, e = def_range(61, "plot_delta_total_and_by_treatment")
fn61 = "\n".join(L(61)[s:e])
fn61 = sub1(fn61, "    savepath=None,", "    save=True,\n    save_path=None,", "c61 sig")
fn61 = sub1(
    fn61,
    """    if savepath:
        fig.savefig(savepath, dpi=dpi, bbox_inches="tight", pad_inches=0.01)
    plt.show()
    return fig, ax
""",
    """    _save(fig, save, save_path)
    return fig, ax
""",
    "c61 tail",
)
fn61 = sub1(fn61, 'cmap = plt.get_cmap("tab10")', 'cmap = mpl.colormaps.get_cmap("tab10")', "c61 cmap")
add(fn61)

# c64 ax helpers (shared by the combined + 3-panel + 2x2 figures)
s, e = def_range(64, "plot_total_elpd_from_summary_ax_black")
blk = "\n".join(L(64)[s:e])
blk = sub1(blk, "def plot_total_elpd_from_summary_ax_black(", "def _ax_total_elpd_black(", "c64 black rename")
add(blk)
s, e = def_range(64, "plot_delta_total_and_by_treatment_ax_colorblind_bracket")
blk = "\n".join(L(64)[s:e])
blk = sub1(blk, "def plot_delta_total_and_by_treatment_ax_colorblind_bracket(", "def _ax_delta_bracket(", "c64 bracket rename")
add(blk)

# c64 combined figure
c64_body = body(64, 172, 215)
c64_body = sub1(c64_body, "plot_total_elpd_from_summary_ax_black(", "_ax_total_elpd_black(", "c64 call a")
c64_body = sub1(c64_body, "plot_delta_total_and_by_treatment_ax_colorblind_bracket(", "_ax_delta_bracket(", "c64 call b")
c64_body = sub1(
    c64_body,
    """savepath = "../figures/bats/psisloo_combined.svg"
os.makedirs(os.path.dirname(savepath), exist_ok=True)
fig.savefig(savepath, bbox_inches="tight", pad_inches=0.01)

plt.show()
""",
    """_save(fig, save, save_path)
return fig
""",
    "c64 save",
)
add(
    "def plot_psisloo_combined(res_total, delta_df, treatment_order, *, save=True, save_path=None):\n"
    "    \"\"\"Two stacked panels: total PSIS-LOO ELPD of both models, and the vMU - vM\n"
    "    ELPD difference total + within-treatment bracket.\"\"\"\n"
    + indent(c64_body)
)

# c66 -> plot_delta_subplots_by_treatment
s, e = def_range(66, "plot_delta_subplots_by_treatment")
fn66 = "\n".join(L(66)[s:e])
fn66 = sub1(
    fn66,
    "def plot_delta_subplots_by_treatment(res_ind, treatment_order, *, savepath=None):",
    "def plot_delta_subplots_by_treatment(res_ind, treatment_order, *, save=True, save_path=None):",
    "c66 sig",
)
fn66 = sub1(
    fn66,
    """    if savepath:
        # save original requested file
        fig.savefig(savepath, dpi=300, bbox_inches="tight")

        # ALSO save SVG with same basename/style
        root, _ = os.path.splitext(savepath)
        fig.savefig(f"{root}.svg", bbox_inches="tight", pad_inches=0.01)

    plt.show()
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c66 tail",
)
add(fn66)

# ============================================================ c71 population mean-HDI figure
s, e = def_range(71, "plot_mean_hdi_horizontal")
fn71 = "\n".join(L(71)[s:e])
fn71 = sub1(fn71, "    xlim=None, savepath=None,", "    xlim=None, save=True, save_path=None,", "c71 sig")
fn71 = drop_line(fn71, "save_svg_also=True", "c71 save_svg_also")
fn71 = sub1(
    fn71,
    """    plt.tight_layout()
    if savepath:
        fig.savefig(savepath, bbox_inches="tight", pad_inches=0.01)
        if save_svg_also:
            root, _ = os.path.splitext(savepath)
            fig.savefig(root + ".svg", bbox_inches="tight", pad_inches=0.01)
    plt.show()
    return fig, ax
""",
    """    plt.tight_layout()
    _save(fig, save, save_path)
    return fig, ax
""",
    "c71 tail",
)
add(fn71)

# ============================================================ c73 2x2 params figure
s, e = def_range(73, "plot_mean_hdi_horizontal_on_ax")
add("\n".join(L(73)[s:e]))

c73_body = body(73, 435, 556)
c73_body = sub1(c73_body, "res_params = summarize_pop_params_for_2x2_by_treatment(", "res_params = posthoc.summarize_pop_params_for_2x2_by_treatment(", "c73 posthoc call")
c73_body = sub1(c73_body, "plot_delta_total_and_by_treatment_ax_colorblind_bracket(", "_ax_delta_bracket(", "c73 bracket call")
c73_body = sub1(
    c73_body,
    """# --- 6) Save ---
savepath = "../figures/bats/psisloo_plus_params_2x2.svg"
os.makedirs(os.path.dirname(savepath), exist_ok=True)
fig.savefig(savepath, bbox_inches="tight", pad_inches=0.01)

# also PDF
root, _ = os.path.splitext(savepath)
fig.savefig(root + ".pdf", bbox_inches="tight", pad_inches=0.01)

plt.show()
""",
    """_save(fig, save, save_path)
return fig
""",
    "c73 save",
)
add(
    "def plot_psisloo_plus_params_2x2(\n"
    "    delta_df, treatment_order, idata_mix_full, df_track,\n"
    "    *, hdi_prob=0.95, max_draws=2000, save=True, save_path=None,\n"
    "):\n"
    "    \"\"\"ELPD-difference panel on top of the full-model population summaries\n"
    "    (mean angle, mixture weight, MRL, directed-component density at 0) in a 2x2 grid.\"\"\"\n"
    + indent(c73_body)
)

# ============================================================ c75 matched predictive density
s, e = def_range(75, "plot_predictive_density_hdi_ax_colorblind")
blk = "\n".join(L(75)[s:e])
blk = sub1(blk, "def plot_predictive_density_hdi_ax_colorblind(", "def _ax_pred_density(", "c75 ax rename")
add(blk)

c75_body = body(75, 262, 377)
c75_body = sub1(
    c75_body,
    "res_treat = build_res_treat_from_vonmises_idata(\n    _idata_vm,",
    "res_treat = posthoc.build_res_treat_from_vonmises_idata(\n    idata_vmn,",
    "c75 posthoc call",
)
c75_body = subn(c75_body, "plot_predictive_density_hdi_ax_colorblind(", "_ax_pred_density(", 2, "c75 ax calls")
c75_body = sub1(
    c75_body,
    """plt.tight_layout()

# save both SVG + PDF (same basename)
savebase = "../figures/bats/predictive_density_two_model_comparison_matched_vmu_directed_vs_vm"
os.makedirs(os.path.dirname(savebase), exist_ok=True)
fig.savefig(savebase + ".svg", bbox_inches="tight", pad_inches=0.01)
fig.savefig(savebase + ".pdf", bbox_inches="tight", pad_inches=0.01)

plt.show()
""",
    """plt.tight_layout()
_save(fig, save, save_path)
return fig
""",
    "c75 save",
)
add(
    "def plot_predictive_density_vmu_directed_vs_vm(res_full, df_track, idata_vmn, *, save=True, save_path=None):\n"
    "    \"\"\"Matched two-panel comparison at the population level: the vMU model's\n"
    "   Directed-component-only predictive density at 0 vs the von Mises model's.\"\"\"\n"
    + indent(c75_body)
)

# ============================================================ c76 three-panel figure
c76_body = body(76, 31, 132)
c76_body = sub1(c76_body, "plot_total_elpd_from_summary_ax_black(", "_ax_total_elpd_black(", "c76 call a")
c76_body = sub1(c76_body, "plot_delta_total_and_by_treatment_ax_colorblind_bracket(", "_ax_delta_bracket(", "c76 call b")
c76_body = sub1(c76_body, "plot_predictive_density_hdi_ax_colorblind(", "_ax_pred_density(", "c76 call c")
c76_body = sub1(
    c76_body,
    """# save SVG + PDF
savebase = "../figures/bats/psisloo_plus_predictive_density_3panel"
os.makedirs(os.path.dirname(savebase), exist_ok=True)
fig.savefig(savebase + ".svg", bbox_inches="tight", pad_inches=0.01)
fig.savefig(savebase + ".pdf", bbox_inches="tight", pad_inches=0.01)

plt.show()
""",
    """_save(fig, save, save_path)
return fig
""",
    "c76 save",
)
add(
    "def plot_psisloo_plus_predictive_density_3panel(res_total, delta_df, treatment_order, res_full, *, save=True, save_path=None):\n"
    "    \"\"\"Three stacked panels: total ELPD, vMU - vM ELPD difference, and the vMU\n"
    "    mixture predictive density at 0 by treatment.\"\"\"\n"
    + indent(c76_body)
)

# ============================================================ c79 tag-level subplots
s, e = def_range(79, "plot_individual_subplots_by_treatment")
fn79 = "\n".join(L(79)[s:e])
fn79 = sub1(fn79, "    savepath=None,", "    save=True,\n    save_path=None,", "c79 sig")
fn79 = drop_line(fn79, "save_svg_also=True", "c79 save_svg_also")
fn79 = sub1(
    fn79,
    """    if savepath:
        fig.savefig(savepath, dpi=DPI, bbox_inches="tight", pad_inches=0.01)
        if save_svg_also:
            root, _ = os.path.splitext(savepath)
            fig.savefig(root + ".svg", bbox_inches="tight", pad_inches=0.01)

    plt.show()
    return fig, axes
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c79 tail",
)
add(fn79)

# ============================================================ c80 population PPC
s80, e80 = topdef(80, "plot_ppc_by_treatment")
c80_body = "\n".join(L(80)[e80 + 1 : 645])  # body lines 518..645 (1-based)
c80_body = sub1(
    c80_body,
    """    fig.savefig(outpath_pdf, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(outpath_svg, bbox_inches="tight", pad_inches=0.02)
    plt.show()
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c80 save",
)
add(
    "def plot_ppc_by_treatment(\n"
    "    treat_levels, angles, bin_centers, obs_boot, ppc, treatment_label_map=None,\n"
    "    *, save=True, save_path=None,\n"
    "):\n"
    "    \"\"\"Polar population-level PPC per treatment: bootstrap KDE bands of the\n"
    "    observed goal-relative headings vs the posterior-predictive KDE bands.\"\"\"\n"
    "    EDGE_RADIUS = 0.55\n"
    "    BAND_ALPHA = 0.25\n"
    "    HDI_PROB = 0.95\n"
    + c80_body
)

# ============================================================ c84 individual PPC grid
c84_body = body(84, 8, 115)
c84_body = sub1(c84_body, "def format_polar(ax, rmax):", "def _polar_format_grid(ax, rmax):", "c84 polar rename")
c84_body = sub1(c84_body, "format_polar(ax, rmax=RMAX)", "_polar_format_grid(ax, rmax=RMAX)", "c84 polar call")
c84_body = sub1(
    c84_body,
    """out_pdf = "../figures/bats/bat_individual_ppc_bootstrap_vs_ppc_kde.pdf"
out_svg = "../figures/bats/bat_individual_ppc_bootstrap_vs_ppc_kde.svg"
fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.005)
fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.005)
plt.show()
""",
    """_save(fig, save, save_path)
return fig
""",
    "c84 save",
)
add(
    "def plot_individual_ppc_grid(ppc, *, save=True, save_path=None):\n"
    "    \"\"\"Polar grid (one panel per tag) of each animal's bootstrap KDE bands vs\n"
    "    its posterior-predictive KDE bands; ``ppc`` is the dict from\n"
    "    ``posthoc.individual_ppc``.\"\"\"\n"
    "    results = ppc[\"results\"]\n"
    "    angles = ppc[\"angles\"]\n"
    "    to_r = ppc[\"to_r\"]\n"
    "    RMAX = ppc[\"RMAX\"]\n"
    "    EDGE_RADIUS = ppc[\"EDGE_RADIUS\"]\n"
    "    AMP = ppc[\"AMP\"]\n"
    "    treat_levels = ppc[\"treat_levels\"]\n"
    "    treat_to_color = ppc[\"treat_to_color\"]\n"
    "    nice_treat_label = ppc[\"nice_treat_label\"]\n"
    + indent(c84_body)
)

# ============================================================ c69 2x2 mixture polar
s, e = def_range(69, "format_polar")
blk = "\n".join(L(69)[s:e])
blk = sub1(blk, "def format_polar(ax, *, rmax):", "def _polar_format_mix(ax, *, rmax):", "c69 polar rename")
add(blk)
s, e = def_range(69, "vonmises_pdf")
blk = "\n".join(L(69)[s:e])
blk = sub1(blk, "def vonmises_pdf(theta, mu, kappa):", "def _vonmises_pdf_mix(theta, mu, kappa):", "c69 vm rename")
add(blk)

c69_body = body(69, 145, 335)
c69_body = re.sub(r"(?<![\w.])vonmises_pdf\(", "_vonmises_pdf_mix(", c69_body)
assert "_vonmises_pdf_mix_mix" not in c69_body
assert c69_body.count("_vonmises_pdf_mix(") == 2, c69_body.count("_vonmises_pdf_mix(")
c69_body = sub1(c69_body, "format_polar(ax, rmax=RMAX)", "_polar_format_mix(ax, rmax=RMAX)", "c69 polar call")
c69_body = sub1(
    c69_body,
    """fig.savefig("../figures/bats/vmum_mixture_full_polar_pdfs_population_only.pdf",
            bbox_inches="tight", pad_inches=0.01)
fig.savefig("../figures/bats/vmum_mixture_full_polar_pdfs_population_only.svg",
            bbox_inches="tight", pad_inches=0.01)
plt.show()
""",
    """_save(fig, save, save_path)
return fig
""",
    "c69 save",
)
add(
    "def plot_mixture_polar_pdfs(idata_mix_full, df_track, *, save=True, save_path=None):\n"
    "    \"\"\"2x2 polar panels of the full vMU mixture's predictive density on the\n"
    "    goal-relative circle, population level, one panel per treatment band.\"\"\"\n"
    "    HDI_PROB = 0.95\n"
    "    MAX_DRAWS = 4000\n"
    "    EDGE_RADIUS = 1.0\n"
    "    AMP = 1.0\n"
    "    BAND_ALPHA = 0.35\n"
    "    KAPPA_UNIFORM = 1e-6\n"
    "    THETA_N = 900\n"
    "    SHOW_LEGEND = False\n"
    + indent(c69_body)
)

# ============================================================ c87 viridis tracks
for name in ("_pick_order_col", "_ensure_step_index", "_make_points_gdf"):
    s, e = def_range(87, name)
    add("\n".join(L(87)[s:e]))
for name, new in (
    ("hide_axes", "hide_axes_pgoal"),
    ("add_scale_bar_dual", "add_scale_bar_dual_pgoal"),
    ("add_compass", "add_compass_pgoal"),
    ("_nice_scale_length_m", "_nice_scale_length_m_pgoal"),
):
    s, e = def_range(87, name)
    blk = "\n".join(L(87)[s:e])
    blk = sub1(blk, f"def {name}(", f"def {new}(", f"c87 {name} rename")
    add(blk)

s, e = def_range(87, "_add_scale_and_compass_bottom_right_equal_pad")
blk = "\n".join(L(87)[s:e])
blk = sub1(blk, "    scale_len_m = _nice_scale_length_m(abs(x1 - x0))", "    scale_len_m = _nice_scale_length_m_pgoal(abs(x1 - x0))", "c87 pad nice")
blk = sub1(blk, "        return add_scale_bar_dual(", "        return add_scale_bar_dual_pgoal(", "c87 pad bar")
blk = sub1(blk, "    add_compass(", "    add_compass_pgoal(", "c87 pad compass")
add(blk)

s, e = def_range(87, "plot_treatment_panels_tracks_colored_by_pgoal")
fn87 = "\n".join(L(87)[s:e])
fn87 = sub1(fn87, "def plot_treatment_panels_tracks_colored_by_pgoal(", "def plot_tracks_pgoal_viridis(", "c87 rename")
fn87 = sub1(
    fn87,
    '    out_pdf="../figures/bats/treatment_panels_tracks_pgoal_viridis.pdf",\n    out_svg="../figures/bats/treatment_panels_tracks_pgoal_viridis.svg",\n',
    "",
    "c87 drop out paths",
)
fn87 = sub1(
    fn87,
    "    cbar_labelsize=9,\n):",
    "    cbar_labelsize=9,\n    save: bool = True,\n    save_path=None,\n):",
    "c87 save kwargs",
)
fn87 = sub1(fn87, "hide_axes(", "hide_axes_pgoal(", "c87 hide call")
fn87 = sub1(
    fn87,
    """    fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.01, dpi=300)
    fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.01, dpi=300)
    plt.show()

    return fig, axes
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c87 tail",
)
add(fn87)

# ============================================================ c89 small bigfont tracks
s, e = def_range(89, "_add_scale_and_compass_bottom_right_equal_pad")
blk = "\n".join(L(89)[s:e])
blk = sub1(blk, "def _add_scale_and_compass_bottom_right_equal_pad(ax, corner_pad_px=14):", "def _add_scale_and_compass_small_pad(ax, corner_pad_px=14):", "c89 pad rename")
blk = sub1(blk, "    scale_len_m = _nice_scale_length_m(abs(x1 - x0))", "    scale_len_m = _nice_scale_length_m_pgoal(abs(x1 - x0))", "c89 pad nice")
blk = sub1(blk, "        return add_scale_bar_dual(", "        return add_scale_bar_dual_pgoal(", "c89 pad bar")
blk = sub1(blk, "    add_compass(", "    add_compass_pgoal(", "c89 pad compass")
add(blk)

s, e = def_range(89, "plot_treatment_panels_tracks_colored_by_pgoal")
fn89 = "\n".join(L(89)[s:e])
fn89 = sub1(fn89, "def plot_treatment_panels_tracks_colored_by_pgoal(", "def plot_tracks_pgoal_small_bigfont(", "c89 rename")
fn89 = sub1(
    fn89,
    '    out_pdf="../figures/bats/treatment_panels_tracks_pgoal_compact_bottomcbar.pdf",\n    out_svg="../figures/bats/treatment_panels_tracks_pgoal_compact_bottomcbar.svg",\n',
    "",
    "c89 drop out paths",
)
fn89 = sub1(
    fn89,
    """    # selective rasterization controls (NEW)
    rasterize_tracks_pdf=True,     # make PDF lighter/faster
    rasterize_tracks_svg=False,    # keep SVG tracks vector by default
""",
    "",
    "c89 drop rasterize",
)
fn89 = sub1(
    fn89,
    "    svg_raster_dpi=600,            # used only if rasterize_tracks_svg=True\n):",
    "    svg_raster_dpi=600,  # raster dpi of a saved SVG export\n    save: bool = True,\n    save_path=None,\n):",
    "c89 save kwargs",
)
fn89 = sub1(fn89, "hide_axes(", "hide_axes_pgoal(", "c89 hide call")
fn89 = sub1(fn89, "_add_scale_and_compass_bottom_right_equal_pad(", "_add_scale_and_compass_small_pad(", "c89 pad call")
fn89 = sub1(
    fn89,
    """    # save with selective rasterization (tracks only)
    if out_pdf:
        _maybe_mkdir_for_file(out_pdf)
        _set_rasterized(track_artists, rasterize_tracks_pdf)
        fig.savefig(out_pdf, bbox_inches="tight", pad_inches=0.01, dpi=pdf_raster_dpi)

    if out_svg:
        _maybe_mkdir_for_file(out_svg)
        _set_rasterized(track_artists, rasterize_tracks_svg)
        fig.savefig(out_svg, bbox_inches="tight", pad_inches=0.01, dpi=svg_raster_dpi)

    plt.show()
    return fig, axes
""",
    """    _rasterize = str(save_path).endswith(".pdf")
    _set_rasterized(track_artists, _rasterize)
    _save(fig, save, save_path, dpi=pdf_raster_dpi if _rasterize else svg_raster_dpi)
    return fig
""",
    "c89 tail",
)
add(fn89)

# ============================================================ c91 population curves (smooth)
s, e = def_range(91, "plot_population_curve_smooth")
fn91 = "\n".join(L(91)[s:e])
fn91 = sub1(fn91, "    savepath=None,", "    save=True,\n    save_path=None,", "c91 sig")
fn91 = drop_line(fn91, "save_svg_also=True", "c91 save_svg_also")
fn91 = sub1(
    fn91,
    """    if savepath:
        fig.savefig(savepath, dpi=DPI, bbox_inches="tight", pad_inches=0.01)
        if save_svg_also:
            root, _ = os.path.splitext(savepath)
            fig.savefig(root + ".svg", bbox_inches="tight", pad_inches=0.01)

    plt.show()
    return fig, ax
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c91 tail",
)
add(fn91)

# ============================================================ c93 compact swap-in curves
s, e = def_range(93, "plot_population_curve_smooth")
fn93 = "\n".join(L(93)[s:e])
fn93 = sub1(fn93, "def plot_population_curve_smooth(", "def plot_population_curve_smooth_compact(", "c93 rename")
fn93 = sub1(fn93, "    line_lw=LINE_LW_DEFAULT,", "    line_lw=1.2,", "c93 lw")
fn93 = sub1(fn93, "    fig_w=SWAPIN_FIG_W,      # long/skinny (~1x4 polar width)", "    fig_w=3.40,", "c93 w")
fn93 = sub1(fn93, "    fig_h=SWAPIN_FIG_H,      # short (~1x4 polar height)", "    fig_h=1.20,", "c93 h")
fn93 = sub1(fn93, "    savepath=None,", "    save=True,\n    save_path=None,", "c93 sig")
fn93 = drop_line(fn93, "save_svg_also=True", "c93 save_svg_also")
fn93 = sub1(
    fn93,
    """    if savepath:
        fig.savefig(savepath, dpi=DPI, bbox_inches="tight", pad_inches=0.01)
        if save_svg_also:
            root, _ = os.path.splitext(savepath)
            fig.savefig(root + ".svg", bbox_inches="tight", pad_inches=0.01)

    plt.show()
    return fig, ax
""",
    """    _save(fig, save, save_path)
    return fig
""",
    "c93 tail",
)
add(fn93)

# ============================================================ c95 faceted tag curves
s, e = def_range(95, "plot_tag_curves_faceted_4panel")
fn95 = "\n".join(L(95)[s:e])
fn95 = sub1(fn95, "    savepath=None,", "    save=True,\n    save_path=None,", "c95 sig")
fn95 = drop_line(fn95, "save_svg_also=True", "c95 save_svg_also")
fn95 = sub1(
    fn95,
    """        if savepath:
            if len(pages) > 1:
                root, ext = os.path.splitext(savepath)
                out_path = f"{root}_p{p_i}{ext}"
            else:
                out_path = savepath

            fig.savefig(out_path, dpi=DPI, bbox_inches="tight", pad_inches=0.01)

            if save_svg_also:
                root_page, _ = os.path.splitext(out_path)
                fig.savefig(f"{root_page}.svg", bbox_inches="tight", pad_inches=0.01)

        plt.show()
        figs.append(fig)
""",
    """        if save:
            root, ext = os.path.splitext(save_path)
            save_page = f"{root}_p{p_i}{ext}" if len(pages) > 1 else save_path
            _save(fig, save, save_page)

        figs.append(fig)
""",
    "c95 tail",
)
add(fn95)

# ============================================================ descriptive per-tag summaries (old module)
old = open("/tmp/opencode/old_bats_plots.py").read()
old = old[old.index("def _treatment_plot_df(") : old.index("# ===================== PORTED CELLS BELOW")]
old = sub1(old, '    cmap = plt.get_cmap("tab10")', '    cmap = mpl.colormaps.get_cmap("tab10")', "old cmap")
# signatures via explicit replacements
old = sub1(old, "def plot_kappa_hat(kappa_df):", "def plot_kappa_hat(kappa_df, *, save=True, save_path=None):", "old sig 1")
old = sub1(old, "def plot_rbar(kappa_df):", "def plot_rbar(kappa_df, *, save=True, save_path=None):", "old sig 2")
old = sub1(old, "def plot_alignment(align_df):", "def plot_alignment(align_df, *, save=True, save_path=None):", "old sig 3")
for stem in (
    "kappa_hat_by_treatment_median_percentiles",
    "rbar_by_treatment_median_percentiles",
    "alignment0_by_treatment_median_percentiles",
):
    old = sub1(
        old,
        f'    _save(fig, "{stem}", dpi=300, pad_inches=0.01)\n    plt.show()\n',
        "    _save(fig, save, save_path)\n    return fig\n",
        f"old save {stem}",
    )
add(old.strip("\n"))

# ---------------------------------------------------------------- assemble
src = "\n\n".join(OUT)
assert not src.lower().startswith("\n")

# ---------------------------------------------------------------- audit
# Residue tokens must not appear in *code*. Docstrings/comments may name the
# conventions they describe (e.g. the header docstring says "followed by
# plt.show()"), so strip docstrings + comment lines before checking.
audit_src = re.sub(r'"""[\s\S]*?"""', "", src)
audit_src = "\n".join(
    l for l in audit_src.split("\n") if not l.lstrip().startswith("#")
)
for bad in (
    "savepath",
    "out_pdf",
    "out_svg",
    "outpath",
    "out_path",
    "save_svg_also",
    "plt.show()",
    "mpl.cm.get_cmap",
    "plt.get_cmap",
    "_maybe_mkdir_for_file",
    "save_dpi",
    "rasterize_trajectories",
    "rasterize_tracks_pdf",
    "rasterize_tracks_svg",
    "SWAPIN_FIG",
    "LINE_LW_DEFAULT",
):
    assert bad not in audit_src, f"residue: {bad}"
# az.hdi's keyword was renamed hdi_prob -> prob in arviz; the local
# circular/linear hdi helpers and the posthoc API keep the hdi_prob name,
# so only flag hdi_prob passed *into az.hdi itself*
assert not re.search(r"az\.hdi\([^)]*?hdi_prob\s*=", audit_src, re.S), "residue: az.hdi(hdi_prob=...)"
expected = [
    "def hide_axes(",
    "def _nice_scale_length_m(",
    "def _pick_order_cols(",
    "def _plot_gdf_capture_artists(",
    "def _set_rasterized(",
    "def _style_trajectory_artists(",
    "def add_scale_bar_dual(",
    "def add_compass(",
    "def _add_scale_and_compass(",
    "def _ensure_dest_cols(",
    "def _dest_label_from_flags(",
    "def _make_segment_and_gap_lines_from_df(",
    "def plot_destination_panels(",
    "def plot_destination_panels_compact(",
    "def _nice_scale_length_m_flow(",
    "def _format_scale_label(",
    "def _fallback_add_scale_and_compass(",
    "def _safe_normalize(",
    "def _make_lines_from_df(",
    "def lake_outer_polygon_largest(",
    "def polygon_mask_on_grid(",
    "def _adaptive_grid_res(",
    "def solve_laplace_dirichlet(",
    "def harmonic_field_on_extent(",
    "def smoothed_distance_field_on_extent(",
    "def roost_field_on_extent(",
    "def plot_destination_flow_fields(",
    "def plot_per_tag_polar_kde(",
    "def plot_treatment_coefficients_mu_kappa(",
    "def plot_treatment_pred_density_at0_mu_kappa(",
    "def plot_total_elpd(",
    "def plot_delta_total_and_by_treatment(",
    "def _ax_total_elpd_black(",
    "def _ax_delta_bracket(",
    "def plot_psisloo_combined(",
    "def plot_delta_subplots_by_treatment(",
    "def plot_mean_hdi_horizontal(",
    "def plot_mean_hdi_horizontal_on_ax(",
    "def plot_psisloo_plus_params_2x2(",
    "def _ax_pred_density(",
    "def plot_predictive_density_vmu_directed_vs_vm(",
    "def plot_psisloo_plus_predictive_density_3panel(",
    "def plot_individual_subplots_by_treatment(",
    "def plot_ppc_by_treatment(",
    "def plot_individual_ppc_grid(",
    "def _polar_format_mix(",
    "def _vonmises_pdf_mix(",
    "def plot_mixture_polar_pdfs(",
    "def _pick_order_col(",
    "def _ensure_step_index(",
    "def _make_points_gdf(",
    "def hide_axes_pgoal(",
    "def add_scale_bar_dual_pgoal(",
    "def add_compass_pgoal(",
    "def _nice_scale_length_m_pgoal(",
    "def _add_scale_and_compass_bottom_right_equal_pad(",
    "def plot_tracks_pgoal_viridis(",
    "def _add_scale_and_compass_small_pad(",
    "def plot_tracks_pgoal_small_bigfont(",
    "def plot_population_curve_smooth(",
    "def plot_population_curve_smooth_compact(",
    "def plot_tag_curves_faceted_4panel(",
    "def _treatment_plot_df(",
    "def _dots_median_interval(",
    "def plot_kappa_hat(",
    "def plot_rbar(",
    "def plot_alignment(",
]
for name in expected:
    assert src.count(name) == 1, f"{name}: {src.count(name)}"
assert src.count("posthoc.") >= 3
compile(src, "goldshtein_bats_plots.py", "exec")

out_path = "/tmp/opencode/wt-unified/notebooks/goldshtein_bats_plots.py"
open(out_path, "w").write(src + "\n")
print(f"OK: wrote {out_path}, {len(src.splitlines())} lines")
