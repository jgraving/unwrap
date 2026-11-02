"""Post-hoc statistics for notebooks/goldshtein_bats.ipynb.

Everything the thin notebook calls in one line: the full-trajectory
modeling table; the descriptive per-tag circular summaries; the PSIS-LOO
comparison machinery between the von Mises-Uniform (vMU) mixture and the
von Mises baseline (totals, per-treatment and per-individual deltas);
the population- and tag-level parameter summaries; the hierarchical
bootstrap posterior predictive checks (population, and one panel per
tag); and the per-step posterior responsibility of the goal component
with its distance-to-destination kernel curves.

The heavy implementations are the annotated notebook cells moved here
verbatim (one copy per duplicated helper family); the statistical odd
choices of the annotated study (e.g. the logit-space p_goal summary)
are preserved.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import arviz as az
import matplotlib as mpl
import scipy.special as sps
from scipy.special import expit, logit
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, List, Tuple
from scipy.stats import norm
from joblib import Parallel, delayed
from tqdm.auto import tqdm

from unwrap.link import circular_modulo

# ---------------------------------------------------------------------------
# Shared numerics (one copy of the per-cell util families; the cells in the
# annotated notebook redefined them locally before every use).
# ---------------------------------------------------------------------------

#: von Mises concentration proxying a flat uniform circle (1/(2*pi) density).
KAPPA_UNIFORM = 1e-6
#: half-width of the "predictive mass around the destination" statistic, deg.
MASS_HALF_WIDTH_DEG = 45.0

def softplus(x):
    x = np.asarray(x, float)
    return np.where(x > 30, x, np.log1p(np.exp(x)))


def sigmoid(x):
    x = np.asarray(x, float)
    return 1.0 / (1.0 + np.exp(-x))


def _stack_samples_first(da):
    da2 = da.stack(sample=("chain", "draw"))
    other = [d for d in da2.dims if d != "sample"]
    da2 = da2.transpose("sample", *other)
    return da2.values.astype(float)


def hdi_1d(x, prob=0.95):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return np.nan, np.nan
    x = np.sort(x)
    m = int(np.floor(prob * n))
    if m < 1:
        return np.nan, np.nan
    widths = x[m:] - x[: n - m]
    j = int(np.argmin(widths))
    return float(x[j]), float(x[j + m])


def hdi_along_samples(arr, prob=0.95):
    arr = np.asarray(arr, float)
    S, T = arr.shape
    lo = np.empty(T, dtype=float)
    hi = np.empty(T, dtype=float)
    for t in range(T):
        lo[t], hi[t] = hdi_1d(arr[:, t], prob=prob)
    return lo, hi



def hdi_along_samples_2col(arr, prob=0.95):
    """HDI across axis=0 for arr shape (S, N) -> (N,2). NaN-safe."""
    arr = np.asarray(arr, float)
    S, N = arr.shape
    out = np.empty((N, 2), float)
    for i in range(N):
        out[i, 0], out[i, 1] = hdi_1d(arr[:, i], prob=prob)
    return out


def circ_mean(x):
    """Circular mean (radians) for samples x."""
    x = np.asarray(x, float)
    return np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x)))


def wrap_0_2pi(x):
    return np.mod(x, 2.0 * np.pi)


def quantile_band(arr, prob=0.95):
    arr = np.asarray(arr, float)
    qlo = (1.0 - prob) / 2.0
    qhi = 1.0 - qlo
    med = np.nanmedian(arr, axis=0)
    lo  = np.nanquantile(arr, qlo, axis=0)
    hi  = np.nanquantile(arr, qhi, axis=0)
    return med, lo, hi


# --- Bessel ratio A(kappa) = I1/I0 (robust) ---
try:
    from scipy.special import i0 as _i0, i1 as _i1

    def bessel_ratio(kappa):
        kappa = np.asarray(kappa, float)
        return _i1(kappa) / _i0(kappa)
except Exception:
    # fallback approximation (good enough for plotting/summary)
    def bessel_ratio(kappa):
        k = np.asarray(kappa, float)
        out = np.empty_like(k)

        small = k < 1e-3
        mid   = (k >= 1e-3) & (k <= 50.0)
        large = k > 50.0

        out[small] = 0.5 * k[small] - (k[small] ** 3) / 16.0
        km = k[mid]
        out[mid] = km / (2.0 + km)
        kl = k[large]
        out[large] = 1.0 - 1.0/(2.0*kl) - 1.0/(8.0*kl*kl)
        return out


def vonmises_pdf_scalar(theta, mu, kappa):
    # theta, mu can be arrays; uses numpy i0
    theta = np.asarray(theta, float)
    mu = np.asarray(mu, float)
    kappa = np.asarray(kappa, float)
    return np.exp(kappa * np.cos(theta - mu)) / (2 * np.pi * np.i0(kappa))


# --- Optional SciPy CDF for circular interval mass ---
try:
    from scipy.stats import vonmises as _sp_vonmises
    _HAS_SCIPY_VONMISES_CDF = True
except Exception:
    _HAS_SCIPY_VONMISES_CDF = False


def vonmises_interval_mass(theta_lo, theta_hi, mu, kappa, *, grid_n=721):
    """
    P(theta_lo <= Theta <= theta_hi) for Von Mises(mu, kappa).

    Uses SciPy CDF when available; falls back to vectorized numerical integration.
    """
    mu = np.asarray(mu, float).reshape(-1)
    kappa = np.asarray(kappa, float).reshape(-1)

    if _HAS_SCIPY_VONMISES_CDF:
        # SciPy's circular cdf can jump by 1 depending on branch; correct by wrapping.
        c_hi = _sp_vonmises.cdf(theta_hi, kappa, loc=mu)
        c_lo = _sp_vonmises.cdf(theta_lo, kappa, loc=mu)
        mass = np.asarray(c_hi - c_lo, float)

        # Bring to [0, 1] for interval length < 2π (our ±45° interval is π/2).
        mass = np.where(mass < 0.0, mass + 1.0, mass)
        mass = np.where(mass > 1.0, mass - 1.0, mass)
        return np.clip(mass, 0.0, 1.0)

    # Fallback: direct integration of pdf over the interval (vectorized across draws)
    theta = np.linspace(theta_lo, theta_hi, int(grid_n), dtype=float)[:, None]  # (G,1)
    mu2 = mu[None, :]       # (1,S)
    k2  = kappa[None, :]    # (1,S)
    pdf = np.exp(k2 * np.cos(theta - mu2)) / (2 * np.pi * np.i0(k2))
    mass = np.trapz(pdf, theta[:, 0], axis=0)
    return np.clip(np.asarray(mass, float), 0.0, 1.0)


def format_treatment_nicely(name: str) -> str:
    s = str(name).strip()
    if not s:
        return "Control"
    parts = [p.strip() for p in s.split(",") if p.strip()]
    def titleish(x):
        return " ".join([w[:1].upper() + w[1:] if w else "" for w in x.lower().split(" ")])
    parts = [titleish(p) for p in parts]
    return ", ".join(parts)


def _format_treatment_name(name: str) -> str:
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return "Unknown"
    s = str(name).strip().lower()
    parts = [p.strip() for p in s.split(",") if p.strip()]
    if "control" in parts and len(parts) > 1:
        parts = [p for p in parts if p != "control"]
    if not parts:
        parts = ["control"]
    parts = sorted(parts)
    parts = [p[:1].upper() + p[1:] for p in parts]
    return ", ".join(parts)


def _treat_levels_and_labels(df_fit, treat_col="treatment", treat_name_col="treatmentName"):
    if isinstance(df_fit[treat_col].dtype, pd.CategoricalDtype):
        treat_levels = list(df_fit[treat_col].cat.categories)
    else:
        treat_levels = sorted(df_fit[treat_col].dropna().unique().tolist())

    if treat_name_col in df_fit.columns:
        name_map = df_fit.groupby(treat_col, observed=False)[treat_name_col].first().to_dict()
        labels = [_format_treatment_name(name_map.get(t, t)) for t in treat_levels]
    else:
        labels = [_format_treatment_name(t) for t in treat_levels]

    return treat_levels, labels


def log_i0_stable(kappa):
    """Stable-ish log(I0(kappa)) with asymptotic fallback."""
    kappa = np.asarray(kappa, float)
    out = np.empty_like(kappa)
    big = kappa > 700
    small = ~big
    out[small] = np.log(np.i0(kappa[small]))
    out[big] = kappa[big] - 0.5 * np.log(2 * np.pi * kappa[big])
    return out


def log_vonmises_pdf(y, mu, kappa):
    """log f(y|mu,kappa) = kappa*cos(y-mu) - log(2π) - log(I0(kappa))"""
    y = np.asarray(y, float)
    mu = np.asarray(mu, float)
    kappa = np.asarray(kappa, float)
    return kappa * np.cos(y - mu) - np.log(2 * np.pi) - log_i0_stable(kappa)


LOG2PI = np.log(2.0 * np.pi)


def _as_sample_obs(arr_da, n_obs, *, want_components=False):
    """
    Convert posterior var to numpy with:
      - if want_components=True: (S, n_obs, 2)
      - else:                   (S, n_obs)
    Works even if dims are unnamed like *_dim_0, *_dim_1.
    """
    da = _stack_samples_first(arr_da)
    dims = [d for d in da.dims if d != "sample"]
    sizes = [da.sizes[d] for d in dims]

    x = da.values.astype(float)

    if not dims:
        # scalar -> broadcast
        if want_components:
            raise ValueError("Got scalar but want components (obs,2).")
        return np.broadcast_to(x.reshape(-1, 1), (x.shape[0], n_obs))

    if want_components:
        if len(dims) != 2:
            raise ValueError(f"Expected 2 dims (obs,comp) for mixture weights; got dims={da.dims} sizes={sizes}")
        # find obs dim (n_obs) and comp dim (2)
        if n_obs not in sizes or 2 not in sizes:
            raise ValueError(
                f"Can't identify obs/comp dims: need sizes include n_obs={n_obs} and 2; got dims={da.dims} sizes={sizes}"
            )
        obs_axis = sizes.index(n_obs)  # among the 2 dims
        comp_axis = sizes.index(2)

        # current x shape is (S, dim0, dim1)
        # we want (S, n_obs, 2)
        if obs_axis == 0 and comp_axis == 1:
            return x
        elif obs_axis == 1 and comp_axis == 0:
            return np.swapaxes(x, 1, 2)
        else:
            # shouldn't happen
            raise ValueError(f"Unexpected axis mapping dims={da.dims} sizes={sizes}")

    # want (S, n_obs)
    if len(dims) == 1:
        if sizes[0] == n_obs:
            return x
        # scalar-ish per tag etc; we can't safely align
        raise ValueError(f"Expected obs dim size {n_obs} but got {sizes[0]} for dims={da.dims}")

    # if more dims, try to find an obs axis and squeeze the rest
    if n_obs in sizes:
        obs_axis = sizes.index(n_obs)  # among dims
        # move obs axis to position 1 (after S)
        x2 = np.moveaxis(x, 1 + obs_axis, 1)
        # now x2 is (S, n_obs, ...)
        # if trailing dims are all 1, squeeze
        trailing = x2.shape[2:]
        if np.prod(trailing) == 1:
            return x2.reshape(x2.shape[0], x2.shape[1])
        raise ValueError(f"Variable has extra non-singleton dims beyond obs: shape={x2.shape} dims={da.dims}")

    raise ValueError(f"Could not find obs dim size {n_obs} in dims={da.dims} sizes={sizes}")


def _softplus(x):
    x = np.asarray(x, float)
    return np.where(x > 30, x, np.log1p(np.exp(x)))


def _log_i0(kappa):
    """
    Stable log(I0(kappa)) for kappa >= 0 using i0e:
      i0(k) = i0e(k) * exp(k)  => log i0(k) = log i0e(k) + k
    """
    kappa = np.asarray(kappa, float)
    # i0e is defined for all real; VM expects kappa >= 0 anyway
    return np.log(sps.i0e(kappa)) + kappa


# ---------------------------------------------------------------------------
# Modeling table
# ---------------------------------------------------------------------------

def build_modeling_table(track_w_treat: pd.DataFrame) -> pd.DataFrame:
    """Build the full-trajectory modeling table (the step-level rows used
    by every model).

    Keeps the usable rows (tag, destination angle error, treatment -- and
    the along-destination distance where present); converts the degree
    angle error ``flightAngDiffDest_deg`` to the wrapped radial
    observation ``HeadingRadians`` in (-pi, pi]; marks ``tag`` and
    ``treatment`` categorical; and standardizes ``log_dist`` into
    ``z_log_dist`` when the distance column is present.
    """
    angle_col = "flightAngDiffDest_deg"   # degrees
    after_col = "afterLocalization"       # 0/1 (or bool)  (kept, NOT used to filter)
    tag_col   = "tag"
    dist_col  = "dist_to_dest_grad_m"     # meters (already computed in df)

    # -------------------------
    # 1) Preprocess (FULL trajectory)
    # -------------------------
    df_track = track_w_treat.copy()

    # Keep only usable rows
    need = [tag_col, angle_col, "treatment"]
    if dist_col in df_track.columns:
        need.append(dist_col)

    df_track = df_track.dropna(subset=need).copy()

    # Degrees -> radians; wrap to (-pi, pi]
    theta = np.deg2rad(np.mod(pd.to_numeric(df_track[angle_col], errors="coerce"), 360.0))
    df_track["HeadingRadians"] = (theta + np.pi) % (2 * np.pi) - np.pi
    df_track = df_track.dropna(subset=["HeadingRadians"]).copy()

    # Categorical columns
    df_track[tag_col] = df_track[tag_col].astype("category")
    df_track["treatment"] = df_track["treatment"].astype("category")

    # Keep afterLocalization around as int if present (not filtering here)
    if after_col in df_track.columns:
        df_track[after_col] = pd.to_numeric(df_track[after_col], errors="coerce").fillna(0).astype(int)

    # Optional: standardize log distance if available
    if dist_col in df_track.columns:
        df_track["log_dist"] = np.log1p(pd.to_numeric(df_track[dist_col], errors="coerce"))
        log_mean = float(df_track["log_dist"].mean())
        log_sd   = float(df_track["log_dist"].std(ddof=0))
        df_track["z_log_dist"] = (df_track["log_dist"] - log_mean) / (log_sd if log_sd > 0 else 1.0)

        print("log_dist mean:", log_mean, "sd:", log_sd)
    else:
        print(f"Note: {dist_col} not found, skipping log_dist / z_log_dist.")

    print("Full-trajectory rows:", len(df_track))
    print("Full-trajectory tags:", df_track[tag_col].nunique())
    print("Treatments:", df_track["treatment"].cat.categories.tolist())

    return df_track

# ---------------------------------------------------------------------------
# Descriptive circular statistics (per tag, full trajectories)
# ---------------------------------------------------------------------------


def estimate_kappa_from_Rbar(Rbar: float, n: int | None = None) -> float:
    """
    Approximate kappa for von Mises from mean resultant length Rbar in [0,1].
    Optionally applies a common small-sample correction when n is provided.
    """
    if Rbar is None or not np.isfinite(Rbar):
        return np.nan
    Rbar = float(np.clip(Rbar, 0.0, 0.999999999))  # avoid 1 exactly

    if Rbar < 0.53:
        kappa = 2*Rbar + Rbar**3 + (5/6)*Rbar**5
    elif Rbar < 0.85:
        kappa = -0.4 + 1.39*Rbar + 0.43/(1 - Rbar)
    else:
        kappa = 1.0 / (Rbar**3 - 4*Rbar**2 + 3*Rbar)

    # small-sample correction (common in circstats implementations)
    if n is not None and np.isfinite(n) and n > 1:
        n = int(n)
        if kappa < 2:
            kappa = max(kappa - 2/(n*kappa), 0.0) if kappa > 0 else 0.0
        else:
            kappa = (n - 1)**3 * kappa / (n**3 + n)

    return float(kappa)


def mean_resultant_length(theta_rad: np.ndarray) -> float:
    """Mean resultant length Rbar for angles in radians."""
    theta = np.asarray(theta_rad, float)
    theta = theta[np.isfinite(theta)]
    if theta.size == 0:
        return np.nan
    C = np.mean(np.cos(theta))
    S = np.mean(np.sin(theta))
    return float(np.hypot(C, S))


def kappa_by_tag_full_trajectory(
    track_w_treat: pd.DataFrame,
    *,
    angle_col="flightAngDiffDest_deg",   # degrees
    tag_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",
) -> pd.DataFrame:
    """
    Per-tag circular concentration on the FULL trajectory (no phase split).

    Returns columns:
      tag, n, Rbar, kappa_hat, treatment, treatmentName
    """
    cols = [tag_col, angle_col]
    if treat_col in track_w_treat.columns:
        cols.append(treat_col)
    if treat_name_col in track_w_treat.columns:
        cols.append(treat_name_col)

    df = track_w_treat[cols].copy()
    df = df.dropna(subset=[tag_col, angle_col])

    # degrees -> radians, wrapped to [0, 2pi)
    df["_theta"] = np.deg2rad(np.mod(df[angle_col].astype(float), 360.0))

    out_rows = []
    for tag, sub in df.groupby(tag_col, sort=False):
        th = sub["_theta"].to_numpy()
        n = int(th.size)
        Rbar = mean_resultant_length(th)
        kappa_hat = estimate_kappa_from_Rbar(Rbar, n=n)

        row = {
            tag_col: tag,
            "n": n,
            "Rbar": Rbar,
            "kappa_hat": kappa_hat,
        }
        if treat_col in sub.columns:
            row[treat_col] = sub[treat_col].iloc[0]
        if treat_name_col in sub.columns:
            row[treat_name_col] = sub[treat_name_col].iloc[0]

        out_rows.append(row)

    out = pd.DataFrame(out_rows)

    # nice ordering
    if treat_col in out.columns:
        out = out.sort_values([treat_col, tag_col]).reset_index(drop=True)
    else:
        out = out.sort_values([tag_col]).reset_index(drop=True)

    return out


def mean_alignment_to_zero(theta_rad: np.ndarray, *, ref_rad: float = 0.0) -> float:
    """
    Alignment to reference direction (default 0): E[cos(theta - ref)] in [-1, 1].
      1   = perfectly aligned with ref
      0   = no net alignment (e.g. roughly uniform / symmetric)
     -1   = perfectly opposite ref
    """
    theta = np.asarray(theta_rad, float)
    theta = theta[np.isfinite(theta)]
    if theta.size == 0:
        return np.nan
    return float(np.mean(np.cos(theta - ref_rad)))


def alignment_by_tag_full_trajectory(
    track_w_treat: pd.DataFrame,
    *,
    angle_col="flightAngDiffDest_deg",   # degrees (diff from destination), so 0 = aligned
    tag_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",
    ref_deg=0.0,
) -> pd.DataFrame:
    """
    Per-tag alignment to ref_deg (default 0) on the FULL trajectory (no phase split).

    Returns columns:
      tag, n, alignment0, treatment, treatmentName
    """
    cols = [tag_col, angle_col]
    if treat_col in track_w_treat.columns:
        cols.append(treat_col)
    if treat_name_col in track_w_treat.columns:
        cols.append(treat_name_col)

    df = track_w_treat[cols].copy().dropna(subset=[tag_col, angle_col])

    # degrees -> radians; wrap to (-pi, pi] (not strictly required for cos, but keeps things tidy)
    theta = np.deg2rad(np.mod(df[angle_col].astype(float), 360.0))
    df["_theta"] = (theta + np.pi) % (2 * np.pi) - np.pi

    ref_rad = np.deg2rad(ref_deg)

    out_rows = []
    for tag, sub in df.groupby(tag_col, sort=False):
        th = sub["_theta"].to_numpy()
        n = int(th.size)

        row = {
            tag_col: tag,
            "n": n,
            "alignment0": mean_alignment_to_zero(th, ref_rad=ref_rad),
        }
        if treat_col in sub.columns:
            row[treat_col] = sub[treat_col].iloc[0]
        if treat_name_col in sub.columns:
            row[treat_name_col] = sub[treat_name_col].iloc[0]

        out_rows.append(row)

    out = pd.DataFrame(out_rows)

    # nice ordering
    if treat_col in out.columns:
        out = out.sort_values([treat_col, tag_col]).reset_index(drop=True)
    else:
        out = out.sort_values([tag_col]).reset_index(drop=True)

    return out



# The percentiles of P = sigmoid(X) for X ~ Normal(mu, sigma): converts
# normal-scale prior quantiles of the mixture-weight logit to probability
# quantiles (used to read off the p_goal posteriors).
def sigmoid_normal_percentiles(mu, sigma, qs=(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99)):
    """
    Percentiles of P = sigmoid(X), where X ~ Normal(mu, sigma).
    qs are quantiles in [0,1].
    """
    z = norm.ppf(qs)              # standard normal quantiles
    xq = mu + sigma * z           # corresponding X quantiles
    pq = expit(xq)                # transform to P quantiles
    return np.asarray(qs), pq


# ---------------------------------------------------------------------------
# PSIS-LOO model comparison (vMU mixture vs the von Mises baseline)
# ---------------------------------------------------------------------------
def graph_free_rvs(idata, model):
    """The fit result with a posterior group containing all of the model's
    free RVs, i.e. evaluable by the PyMC graph.

    The kappa random effects are fitted non-centered: in the PyMC graph the
    free RVs are then ``*_sigma`` (scale) and ``*_offset`` (factor), but
    bambi stores the transformed centered effects and the ``*_sigma`` scales
    and drops the ``*_offset`` factors. ``pm.compute_deterministics`` and
    ``pm.compute_log_likelihood`` both evaluate the graph against every free
    RV, so the missing factors are reconstructed: the graph computes
    ``effect = sigma * offset`` exactly (sigma scalar per draw, offset one
    factor per group), hence ``offset = effect / sigma`` (0 where sigma is
    0).
    """
    posterior = idata.posterior
    missing = [
        rv.name
        for rv in model.backend.model.free_RVs
        if rv.name not in posterior.data_vars
    ]
    if not missing:
        return idata
    aug = idata.copy(deep=False)
    ds = posterior.to_dataset()
    for name in missing:
        sigma = ds[name.replace("_offset", "_sigma")]
        effect = ds[name[: -len("_offset")]]
        # The graph computes effect = sigma * offset (sigma scalar per
        # draw, offset one factor per group), so the offset carries the
        # group dimension of the centered effect.
        e = np.asarray(effect)
        s = np.asarray(sigma).reshape(np.asarray(sigma).shape + (1,) * (e.ndim - np.asarray(sigma).ndim))
        offset = np.divide(e, s, out=np.zeros(e.shape), where=s != 0)
        ds[name] = (tuple(effect.dims), offset)
    aug["posterior"] = ds
    return aug


def _loo_pointwise_elpd(idata, *, var_name=None):
    loo = az.loo(idata, pointwise=True, var_name=var_name)
    loo_i = np.asarray(loo.loo_i).reshape(-1)
    k = np.asarray(loo.pareto_k).reshape(-1)
    return loo_i, k, loo


def _se_sum(v: np.ndarray) -> float:
    """Rough SE for a SUM: sqrt(n * var(v))."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    n = v.size
    if n <= 1:
        return np.nan
    return float(np.sqrt(n * np.var(v, ddof=1)))


def _pareto_report(k: np.ndarray, threshold: float = 0.7) -> dict:
    k = np.asarray(k, float)
    k = k[np.isfinite(k)]
    if k.size == 0:
        return {"n": 0, "n_gt": 0, "max_k": np.nan, "mean_excess": np.nan, "max_excess": np.nan}

    excess = k - threshold
    gt = excess > 0
    n_gt = int(np.sum(gt))

    if n_gt == 0:
        return {
            "n": int(k.size),
            "n_gt": 0,
            "max_k": float(np.max(k)),
            "mean_excess": 0.0,
            "max_excess": 0.0,
        }

    return {
        "n": int(k.size),
        "n_gt": n_gt,
        "max_k": float(np.max(k)),
        "mean_excess": float(np.mean(excess[gt])),
        "max_excess": float(np.max(excess[gt])),
    }


def compute_psisloo_elpd_summary(
    models,          # dict: {name: idata}  (insertion order matters)
    df_fit,          # EXACT df used in all fits (same rows/order)
    *,
    var_name=None,
    pareto_thresh=0.7,
    order_best_first=True,
):
    """
    Compute PSIS-LOO pointwise ELPD and return a summary table with *absolute* ELPD
    (no deltas).

    Columns:
      - elpd: sum of pointwise elpd (loo_i)
      - se_elpd: rough SE for that sum (based on pointwise variability)
      - is_best: True for max elpd model (optional ordering)
      - pareto diagnostics summary
    """
    if not isinstance(models, dict) or len(models) < 1:
        raise ValueError("models must be a non-empty dict: {name: idata}")

    model_order = list(models.keys())  # preserves insertion order
    n = len(df_fit)

    loo_i_map = {}
    pareto_k_map = {}
    loo_obj_map = {}
    total_elpd = {}
    se_elpd = {}

    for name in model_order:
        idata = models[name]
        loo_i, k, loo = _loo_pointwise_elpd(idata, var_name=var_name)

        if loo_i.size != n:
            raise ValueError(
                f"Row mismatch for model '{name}': df_fit has {n} rows but loo_i has {loo_i.size}. "
                "df_fit MUST be the exact df used in all fits (same rows/order)."
            )

        loo_i_map[name] = loo_i
        pareto_k_map[name] = k
        loo_obj_map[name] = loo

        total_elpd[name] = float(np.sum(loo_i))
        se_elpd[name] = _se_sum(loo_i)

    best_name = max(total_elpd, key=total_elpd.get)

    if order_best_first:
        ordered_names = [best_name] + [m for m in model_order if m != best_name]
    else:
        ordered_names = model_order

    rows = []
    for name in ordered_names:
        rep = _pareto_report(pareto_k_map[name], threshold=pareto_thresh)
        rows.append({
            "model": name,
            "elpd": total_elpd[name],
            "se_elpd": se_elpd[name],
            "is_best": (name == best_name),
            "pareto_n": rep["n"],
            f"pareto_n_gt_{pareto_thresh}": rep["n_gt"],
            "pareto_max_k": rep["max_k"],
            f"pareto_mean_excess_over_{pareto_thresh}": rep["mean_excess"],
            f"pareto_max_excess_over_{pareto_thresh}": rep["max_excess"],
        })

    out = pd.DataFrame(rows)

    cache = {
        "loo_i_map": loo_i_map,
        "pareto_k_map": pareto_k_map,
        "loo_obj_map": loo_obj_map,
        "total_elpd": total_elpd,
        "se_elpd": se_elpd,
        "model_order": model_order,
        "ordered_names": ordered_names,
        "best_name": best_name,
    }
    return out, best_name, cache


def elpd_diff_total_and_by_treatment(
    idata_vmu,
    idata_vm,
    df_fit,
    *,
    treat_col="treatment",
    treat_name_col="treatmentName",
    var_name=None,
):
    """
    Compute ΔELPD = vMU - vM:
      - total across all rows
      - stratified by treatment group
    df_fit MUST be the exact rows/order used in both fits.
    """
    vmu_i, vmu_k = _loo_pointwise_elpd(idata_vmu, var_name=var_name)
    vm_i,  vm_k  = _loo_pointwise_elpd(idata_vm,  var_name=var_name)

    n = len(df_fit)
    if vmu_i.size != n or vm_i.size != n:
        raise ValueError(
            f"Row mismatch: df_fit has {n} rows but loo_i has "
            f"{vmu_i.size} (vMU) and {vm_i.size} (vM)."
        )

    delta_i = vmu_i - vm_i
    k_max_i = np.maximum(vmu_k, vm_k)

    cols = [treat_col]
    if treat_name_col in df_fit.columns:
        cols.append(treat_name_col)

    tmp = df_fit[cols].copy().reset_index(drop=True)
    tmp["delta_i"] = delta_i
    tmp["k_max_i"] = k_max_i

    rows = []

    # TOTAL row
    rows.append({
        "row_type": "total",
        "treatment": None,
        "label": "Total",
        "n_obs": int(len(tmp)),
        "delta_elpd": float(np.sum(tmp["delta_i"].to_numpy())),
        "se_delta_elpd": _se_sum(tmp["delta_i"].to_numpy()),
        "pareto_k_max": float(np.max(tmp["k_max_i"].to_numpy())) if len(tmp) else np.nan,
    })

    # BY TREATMENT rows
    for treat, sub in tmp.groupby(treat_col, sort=False):
        tname = sub[treat_name_col].iloc[0] if treat_name_col in sub.columns else str(treat)
        rows.append({
            "row_type": "treatment",
            "treatment": treat,
            "label": _format_treatment_name(tname),
            "n_obs": int(len(sub)),
            "delta_elpd": float(np.sum(sub["delta_i"].to_numpy())),
            "se_delta_elpd": _se_sum(sub["delta_i"].to_numpy()),
            "pareto_k_max": float(np.max(sub["k_max_i"].to_numpy())) if len(sub) else np.nan,
        })

    return pd.DataFrame(rows)


def _guess_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


def delta_by_individual_treatment(
    idata_vmu,
    idata_vm,
    df_fit,
    *,
    individual_col=None,
    treat_col="treatment",
    treat_name_col="treatmentName",
    var_name=None,
):
    """
    Pointwise Δelpd = vMU - vM, aggregated per (individual, treatment).
    df_fit MUST be the exact rows/order used to fit BOTH models.
    """
    if individual_col is None:
        individual_col = _guess_col(df_fit, ["individual", "participant", "subject", "id", "ID"])
        if individual_col is None:
            raise ValueError(
                "Couldn't guess individual column. Pass individual_col='...'. "
                f"Available columns: {list(df_fit.columns)}"
            )

    vmu_i, vmu_k = _loo_pointwise_elpd(idata_vmu, var_name=var_name)
    vm_i,  vm_k  = _loo_pointwise_elpd(idata_vm,  var_name=var_name)

    n = len(df_fit)
    if vmu_i.size != n or vm_i.size != n:
        raise ValueError(
            f"Row mismatch: df_fit has {n} rows but loo_i has "
            f"{vmu_i.size} (vMU) and {vm_i.size} (vM)."
        )

    delta_i = vmu_i - vm_i
    k_max_i = np.maximum(vmu_k, vm_k)

    cols = [individual_col, treat_col]
    if treat_name_col in df_fit.columns:
        cols.append(treat_name_col)

    tmp = df_fit[cols].copy().reset_index(drop=True)
    tmp["delta_i"] = delta_i
    tmp["k_max_i"] = k_max_i

    out = []
    gb_cols = [treat_col, individual_col]
    for (treat, indiv), sub in tmp.groupby(gb_cols, sort=False):
        v = sub["delta_i"].to_numpy()
        tname = sub[treat_name_col].iloc[0] if treat_name_col in sub.columns else str(treat)
        out.append({
            "treatment": treat,
            "treatment_label": _format_treatment_name(tname),
            "individual": indiv,
            "n_obs": int(v.size),
            "delta_elpd": float(np.sum(v)),
            "se_delta_elpd": _se_sum(v),
            "pareto_k_max": float(np.max(sub["k_max_i"].to_numpy())) if len(sub) else np.nan,
        })

    return pd.DataFrame(out)


# ---------------------------------------------------------------------------
# Population-level model summaries (per-treatment parameter posteriors)
# ---------------------------------------------------------------------------


def summarize_pop_p_kappa_mu_density0_mass45_alignment0_by_treatment(
    idata_mix_full,
    df_fit,
    *,
    treat_col="treatment",
    treat_name_col="treatmentName",
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=KAPPA_UNIFORM,
    mass_half_width_deg=MASS_HALF_WIDTH_DEG,
    mass_grid_n=721,
):
    """
    Summarize treatment-level posterior quantities for the full mixture model.

    Returns a DataFrame with:
      - p_*           : mixture weight (directed-component probability)
      - kappa_*       : directed-component concentration
      - dens0_mix_*   : full mixture predictive density at theta=0
      - dens0_dir_*   : directed-component-only predictive density at theta=0
      - dens0_*       : backward-compatible alias to dens0_mix_*
      - mass45_*      : full mixture predictive mass within ±mass_half_width_deg
      - align0_*      : full mixture E[cos(theta)]
    """
    post = idata_mix_full.posterior

    for v in ["mu0", "logit_p0", "eta_kappa0"]:
        if v not in post.data_vars:
            raise KeyError(f"Missing '{v}'. Available: {list(post.data_vars)[:120]} ...")

    has_beta = "beta_treat" in post.data_vars
    has_kt   = "kappa_treat" in post.data_vars
    has_mt   = "mu_treat" in post.data_vars

    treat_levels, labels = _treat_levels_and_labels(df_fit, treat_col=treat_col, treat_name_col=treat_name_col)
    ref = treat_levels[0]
    nonref = treat_levels[1:]

    rng = np.random.default_rng(123)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]).reshape(-1))   # (S,)
    lp0 = _stack_samples_first(post["logit_p0"]).reshape(-1)               # (S,)
    ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)             # (S,)

    if has_beta:
        beta = _stack_samples_first(post["beta_treat"])
        if beta.ndim == 1:
            beta = beta[:, None]

    if has_kt:
        kt = _stack_samples_first(post["kappa_treat"])
        if kt.ndim == 1:
            kt = kt[:, None]

    if has_mt:
        mt = _stack_samples_first(post["mu_treat"])
        if mt.ndim == 1:
            mt = mt[:, None]

    S = lp0.size
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=max_draws, replace=False)
        mu0 = mu0[idx]
        lp0 = lp0[idx]
        ek0 = ek0[idx]
        if has_beta: beta = beta[idx, :]
        if has_kt:   kt   = kt[idx, :]
        if has_mt:   mt   = mt[idx, :]
        S = max_draws

    theta_half = np.deg2rad(float(mass_half_width_deg))
    theta_lo, theta_hi = -theta_half, theta_half

    f0_unif = float(vonmises_pdf_scalar(0.0, 0.0, kappa_uniform))  # ~ 1/(2π)
    A_unif  = float(bessel_ratio(kappa_uniform))
    mass_unif = float((theta_hi - theta_lo) / (2.0 * np.pi))        # exact uniform mass in interval

    rows = []
    for i, t in enumerate(treat_levels):
        is_ref = (t == ref)

        # Mixture weight p
        lp = lp0.copy()
        if (not is_ref) and has_beta:
            j = nonref.index(t)
            lp = lp + beta[:, j]
        p = sigmoid(lp)

        # Directed concentration kappa
        ek = ek0.copy()
        if (not is_ref) and has_kt:
            j = nonref.index(t)
            ek = ek + kt[:, j]
        kappa = softplus(ek)

        # Directed mean mu
        mu = mu0.copy()
        if (not is_ref) and has_mt:
            j = nonref.index(t)
            mu = circular_modulo(mu + mt[:, j])

        # Predictive density at theta=0
        # Directed component only (conditional on directed state)
        f0_dir = vonmises_pdf_scalar(0.0, mu, kappa)

        # Full mixture predictive density
        f0_mix = p * f0_dir + (1.0 - p) * f0_unif

        # Predictive mass within ±mass_half_width_deg (full mixture)
        mass_goal = vonmises_interval_mass(theta_lo, theta_hi, mu, kappa, grid_n=mass_grid_n)
        mass_pm = p * mass_goal + (1.0 - p) * mass_unif

        # Alignment to 0 (full mixture)
        A_goal = bessel_ratio(kappa)
        align0 = p * (A_goal * np.cos(mu)) + (1.0 - p) * A_unif

        # Summaries
        p_lo, p_hi = hdi_1d(p, prob=hdi_prob)
        k_lo, k_hi = hdi_1d(kappa, prob=hdi_prob)

        f_mix_lo, f_mix_hi = hdi_1d(f0_mix, prob=hdi_prob)
        f_dir_lo, f_dir_hi = hdi_1d(f0_dir, prob=hdi_prob)

        m_lo, m_hi = hdi_1d(mass_pm, prob=hdi_prob)
        a_lo, a_hi = hdi_1d(align0, prob=hdi_prob)

        rows.append({
            "treatment": t,
            "label": labels[i],

            "p_mean": float(np.mean(p)),
            "p_lo": p_lo,
            "p_hi": p_hi,

            "kappa_mean": float(np.mean(kappa)),
            "kappa_lo": k_lo,
            "kappa_hi": k_hi,

            # Full mixture predictive density at theta=0
            "dens0_mix_mean": float(np.mean(f0_mix)),
            "dens0_mix_lo": f_mix_lo,
            "dens0_mix_hi": f_mix_hi,

            # Directed-component-only predictive density at theta=0
            "dens0_dir_mean": float(np.mean(f0_dir)),
            "dens0_dir_lo": f_dir_lo,
            "dens0_dir_hi": f_dir_hi,

            # Backward-compatible alias (old columns = full mixture)
            "dens0_mean": float(np.mean(f0_mix)),
            "dens0_lo": f_mix_lo,
            "dens0_hi": f_mix_hi,

            "mass45_mean": float(np.mean(mass_pm)),
            "mass45_lo": m_lo,
            "mass45_hi": m_hi,

            "align0_mean": float(np.mean(align0)),
            "align0_lo": a_lo,
            "align0_hi": a_hi,
        })

    out = pd.DataFrame(rows)
    out.attrs["uniform_density_at_0"] = float(f0_unif)
    out.attrs["uniform_mass_pm_deg"] = float(mass_unif)
    out.attrs["mass_half_width_deg"] = float(mass_half_width_deg)
    return out


def summarize_pop_params_for_2x2_by_treatment(
    idata_mix_full,
    df_fit,
    *,
    treat_col="treatment",
    treat_name_col="treatmentName",
    treatment_order=None,
    hdi_prob=0.95,
    max_draws=2000,
):
    """
    Returns DataFrame with columns:
      mu_deg_mean, mu_deg_lo, mu_deg_hi  (mean angle, degrees; HDI computed around circular mean)
      p_mean, p_lo, p_hi                 (mixture weight)
      mrl_mean, mrl_lo, mrl_hi           (I1/I0(kappa))
      dens0_dir_mean, dens0_dir_lo, dens0_dir_hi  (directed-component density at 0)
    """
    post = idata_mix_full.posterior
    for v in ["mu0", "logit_p0", "eta_kappa0"]:
        if v not in post.data_vars:
            raise KeyError(f"Missing '{v}' in posterior. Found: {list(post.data_vars)[:80]} ...")

    has_beta = "beta_treat" in post.data_vars
    has_kt   = "kappa_treat" in post.data_vars
    has_mt   = "mu_treat" in post.data_vars

    treat_levels, labels = _treat_levels_and_labels(df_fit, treat_col=treat_col, treat_name_col=treat_name_col)

    # enforce explicit order if provided
    if treatment_order is not None:
        treat_levels = list(treatment_order)
        # map labels accordingly (fallback to formatted level)
        lab_map = dict(zip(*_treat_levels_and_labels(df_fit, treat_col=treat_col, treat_name_col=treat_name_col)))
        labels = [_format_treatment_name(lab_map.get(t, t)) for t in treat_levels]

    ref = treat_levels[0]
    nonref = treat_levels[1:]

    rng = np.random.default_rng(123)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]).reshape(-1))
    lp0 = _stack_samples_first(post["logit_p0"]).reshape(-1)
    ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)

    if has_beta:
        beta = _stack_samples_first(post["beta_treat"])
        if beta.ndim == 1:
            beta = beta[:, None]
    if has_kt:
        kt = _stack_samples_first(post["kappa_treat"])
        if kt.ndim == 1:
            kt = kt[:, None]
    if has_mt:
        mt = _stack_samples_first(post["mu_treat"])
        if mt.ndim == 1:
            mt = mt[:, None]

    S = lp0.size
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=max_draws, replace=False)
        mu0 = mu0[idx]
        lp0 = lp0[idx]
        ek0 = ek0[idx]
        if has_beta: beta = beta[idx, :]
        if has_kt:   kt   = kt[idx, :]
        if has_mt:   mt   = mt[idx, :]
        S = max_draws

    rows = []
    for i, t in enumerate(treat_levels):
        is_ref = (t == ref)

        # mixture weight
        lp = lp0.copy()
        if (not is_ref) and has_beta:
            j = nonref.index(t)
            lp = lp + beta[:, j]
        p = sigmoid(lp)

        # kappa and MRL
        ek = ek0.copy()
        if (not is_ref) and has_kt:
            j = nonref.index(t)
            ek = ek + kt[:, j]
        kappa = softplus(ek)
        mrl = bessel_ratio(kappa)

        # mu (radians)
        mu = mu0.copy()
        if (not is_ref) and has_mt:
            j = nonref.index(t)
            mu = circular_modulo(mu + mt[:, j])

        # circular mean + HDI around the mean (avoid wrap issues)
        mu_bar = circ_mean(mu)
        dmu = circular_modulo(mu - mu_bar)            # centered deviations in (-pi, pi]
        dmu_deg = np.rad2deg(dmu)
        dlo, dhi = hdi_1d(dmu_deg, prob=hdi_prob)

        mu_bar_deg = float(np.rad2deg(mu_bar))
        mu_lo_deg = mu_bar_deg + float(dlo)
        mu_hi_deg = mu_bar_deg + float(dhi)

        # directed-component predictive density at 0
        f0_dir = vonmises_pdf_scalar(0.0, mu, kappa)

        # summarize
        p_lo, p_hi = hdi_1d(p, prob=hdi_prob)
        mrl_lo, mrl_hi = hdi_1d(mrl, prob=hdi_prob)
        f_lo, f_hi = hdi_1d(f0_dir, prob=hdi_prob)

        rows.append({
            "treatment": t,
            "label": labels[i],

            "mu_deg_mean": mu_bar_deg,
            "mu_deg_lo": mu_lo_deg,
            "mu_deg_hi": mu_hi_deg,

            "p_mean": float(np.mean(p)),
            "p_lo": p_lo,
            "p_hi": p_hi,

            "mrl_mean": float(np.mean(mrl)),
            "mrl_lo": mrl_lo,
            "mrl_hi": mrl_hi,

            "dens0_dir_mean": float(np.mean(f0_dir)),
            "dens0_dir_lo": f_lo,
            "dens0_dir_hi": f_hi,
        })

    return pd.DataFrame(rows)


def build_res_treat_from_vonmises_idata(
    idata_treatment,
    df_track,
    *,
    theta0_deg=0.0,
    hdi_prob=0.95,
):
    post = idata_treatment.posterior
    theta0 = np.deg2rad(theta0_deg)

    # ---------- helpers ----------
    def _has(varname: str) -> bool:
        return varname in post.data_vars

    def circular_modulo(x):
        return (x + np.pi) % (2 * np.pi) - np.pi

    def softplus(x):
        x = np.asarray(x, float)
        return np.where(x > 30, x, np.log1p(np.exp(x)))

    def vonmises_pdf(theta, mu, kappa):
        # pdf(theta | mu, kappa) = exp(kappa cos(theta-mu)) / (2π I0(kappa))
        return np.exp(kappa * np.cos(theta - mu)) / (2.0 * np.pi * np.i0(kappa))

    def _stack_samples(da):
        return da.stack(sample=("chain", "draw")).values.astype(float)

    def _format_treatment_name(name: str):
        s = str(name).strip().lower()
        parts = [p.strip() for p in s.split(",") if p.strip()]
        if "control" in parts and len(parts) > 1:
            parts = [p for p in parts if p != "control"]
        if not parts:
            parts = ["control"]
        parts = sorted(parts)
        parts = [p[:1].upper() + p[1:] for p in parts]
        return ", ".join(parts)

    def _cat_dim(varname: str):
        da = post[varname]
        dims = [d for d in da.dims if d not in ("chain", "draw")]
        if len(dims) != 1:
            raise ValueError(f"Expected 1 non-(chain/draw) dim for '{varname}', got {da.dims}")
        return dims[0]

    def _sel_level(varname: str, lvl, nonref_levels):
        dim = _cat_dim(varname)
        coords = list(post[varname].coords[dim].values)
        if lvl in coords:
            return {dim: lvl}
        idx = nonref_levels.index(lvl)
        return {dim: coords[idx]}

    # ---------- require intercepts ----------
    if not _has("Intercept"):
        raise KeyError(f"Missing 'Intercept' (needed for μ). Available: {list(post.data_vars)[:40]} ...")
    if not _has("kappa_Intercept"):
        raise KeyError(f"Missing 'kappa_Intercept' (needed for κ). Available: {list(post.data_vars)[:40]} ...")

    mu_has_treat = _has("treatment")
    k_has_treat  = _has("kappa_treatment")

    # ---------- treatment levels ----------
    treat_levels = (
        df_track["treatment"].cat.categories.tolist()
        if hasattr(df_track["treatment"], "cat")
        else sorted(df_track["treatment"].unique())
    )
    ref = treat_levels[0]
    nonref = treat_levels[1:]

    # ---------- pretty labels ----------
    nice = {t: str(t) for t in treat_levels}
    if "treatmentName" in df_track.columns:
        name_map = df_track.groupby("treatment", observed=False)["treatmentName"].first().to_dict()
        for t in treat_levels:
            nice[t] = _format_treatment_name(name_map.get(t, t))
    else:
        for t in treat_levels:
            nice[t] = _format_treatment_name(t)

    # ---------- rows to plot ----------
    rows_have_treat = (mu_has_treat or k_has_treat) and (len(nonref) > 0)
    if rows_have_treat:
        labels = ["Control"] + [nice[t] for t in nonref]
        row_treatment = [ref] + nonref
    else:
        labels = ["Control"]
        row_treatment = [ref]

    # ---------- population-level intercept draws ----------
    b0_mu = _stack_samples(post["Intercept"])        # (S,)
    b0_k  = _stack_samples(post["kappa_Intercept"])  # (S,)

    # ---------- compute predictive density at theta0 ----------
    dens_draws_by_row = []
    for t in row_treatment:
        eta_mu = b0_mu.copy()
        eta_k  = b0_k.copy()

        if t != ref:
            if mu_has_treat:
                b_trt_mu = _stack_samples(post["treatment"].sel(**_sel_level("treatment", t, nonref)))
                eta_mu = eta_mu + b_trt_mu
            if k_has_treat:
                b_trt_k = _stack_samples(post["kappa_treatment"].sel(**_sel_level("kappa_treatment", t, nonref)))
                eta_k = eta_k + b_trt_k

        mu    = circular_modulo(eta_mu)
        kappa = softplus(eta_k)
        dens0 = vonmises_pdf(theta0, mu, kappa)  # (S,)
        dens_draws_by_row.append(dens0)

    # ---------- summaries ----------
    dens0_mean = np.array([d.mean() for d in dens_draws_by_row], float)
    hdi = np.array([az.hdi(d, prob=hdi_prob) for d in dens_draws_by_row], float)  # (K, 2)
    dens0_lo = hdi[:, 0]
    dens0_hi = hdi[:, 1]

    res_treat = pd.DataFrame({
        "label": labels,
        "dens0_mean": dens0_mean,
        "dens0_lo": dens0_lo,
        "dens0_hi": dens0_hi,
    })

    return res_treat


def summarize_taglevel_from_model_terms(
    idata_mix_full,
    df_fit,
    *,
    treat_col="treatment",
    treat_name_col="treatmentName",
    individual_col="tag",          # your "individual" is tag
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=KAPPA_UNIFORM,
    seed=123,
):
    """
    Build tag-level params from pop + treatment + tag offsets:
      mu_it    = wrap(mu0 + mu_treat[t] + mu_tag[i])
      kappa_it = softplus(eta_kappa0 + kappa_treat[t] + eta_kappa_tag[i] * sigma_kappa_tag)
      p_it     = sigmoid(logit_p0 + beta_treat[t] + logit_p_tag[i] * sigma_p_tag)

    Then derived:
      dens0_it, align0_it
    """
    post = idata_mix_full.posterior

    # required vars in *your* model
    required = ["mu0", "mu_tag", "eta_kappa0", "eta_kappa_tag", "sigma_kappa_tag",
                "logit_p0", "logit_p_tag", "sigma_p_tag"]
    missing = [v for v in required if v not in post.data_vars]
    if missing:
        raise KeyError(f"Missing required posterior vars: {missing}. Available: {list(post.data_vars)[:120]} ...")

    # optional treatment contrasts
    has_mt   = "mu_treat" in post.data_vars
    has_kt   = "kappa_treat" in post.data_vars
    has_beta = "beta_treat" in post.data_vars

    # treatment levels
    treat_levels, treat_labels = _treat_levels_and_labels(df_fit, treat_col=treat_col, treat_name_col=treat_name_col)
    n_treat = len(treat_levels)

    ref = treat_levels[0]
    nonref = treat_levels[1:]
    treat_to_j = {t: nonref.index(t) for t in nonref}  # nonref -> 0..K-1

    # map each tag -> its (single) treatment (mode if messy)
    tmp = df_fit[[individual_col, treat_col]].copy()
    tmp[treat_col] = pd.Categorical(tmp[treat_col], categories=treat_levels, ordered=True)
    tag_to_treat = (tmp.groupby(individual_col, observed=False)[treat_col]
                      .agg(lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0])
                      .to_dict())

    tag_codes = np.array(sorted(tag_to_treat.keys(), key=lambda x: str(x)), dtype=object)
    n_tag = len(tag_codes)

    # IMPORTANT: tag index order must match the model's tag_idx ordering.
    # If you built tag_idx with np.unique(df_track["tag"], return_inverse=True), then
    # tag_codes_model is np.unique(df_track["tag"])[0] in sorted order.
    # We'll align to that by using np.unique on df_fit (same operation).
    tag_codes_model = np.unique(df_fit[individual_col].to_numpy())
    if len(tag_codes_model) != n_tag:
        raise ValueError("Tag code mismatch between df grouping and np.unique. Check df contents.")

    # posterior arrays
    rng = np.random.default_rng(seed)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]).reshape(-1))                 # (S,)
    mu_tag = _stack_samples_first(post["mu_tag"])                                        # (S, n_tag)

    ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)                           # (S,)
    ek_tag = _stack_samples_first(post["eta_kappa_tag"])                                 # (S, n_tag)
    sigk = _stack_samples_first(post["sigma_kappa_tag"]).reshape(-1)                     # (S,)

    lp0 = _stack_samples_first(post["logit_p0"]).reshape(-1)                             # (S,)
    lp_tag = _stack_samples_first(post["logit_p_tag"])                                   # (S, n_tag)
    sigp = _stack_samples_first(post["sigma_p_tag"]).reshape(-1)                         # (S,)

    # treatment arrays -> full (S, n_treat) with 0 at ref
    def _full_contr(arr_name):
        if arr_name not in post.data_vars:
            return None
        arr = _stack_samples_first(post[arr_name])
        if arr.ndim == 1:
            arr = arr[:, None]
        if arr.shape[1] != (n_treat - 1):
            raise ValueError(f"{arr_name} second dim {arr.shape[1]} != n_treat-1={n_treat-1}")
        z = np.zeros((arr.shape[0], 1), float)
        return np.concatenate([z, arr], axis=1)  # (S, n_treat)

    mu_treat_full = _full_contr("mu_treat")       # (S, n_treat) or None
    ktreat_full   = _full_contr("kappa_treat")    # (S, n_treat) or None
    beta_full     = _full_contr("beta_treat")     # (S, n_treat) or None

    # basic shape checks
    S = mu0.size
    if mu_tag.shape[1] != len(tag_codes_model):
        raise ValueError(f"mu_tag has n={mu_tag.shape[1]} but df unique tags has {len(tag_codes_model)}.")
    if ek_tag.shape[1] != len(tag_codes_model) or lp_tag.shape[1] != len(tag_codes_model):
        raise ValueError("Tag RE shapes don't match df unique tags; check tag indexing order.")

    # subsample draws once
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=max_draws, replace=False)
        mu0 = mu0[idx]
        mu_tag = mu_tag[idx, :]
        ek0 = ek0[idx]
        ek_tag = ek_tag[idx, :]
        sigk = sigk[idx]
        lp0 = lp0[idx]
        lp_tag = lp_tag[idx, :]
        sigp = sigp[idx]
        if mu_treat_full is not None: mu_treat_full = mu_treat_full[idx, :]
        if ktreat_full   is not None: ktreat_full   = ktreat_full[idx, :]
        if beta_full     is not None: beta_full     = beta_full[idx, :]
        S = max_draws

    # constants for uniform component
    f0_unif = float(vonmises_pdf_scalar(0.0, 0.0, kappa_uniform))
    A_unif  = float(bessel_ratio(kappa_uniform))

    # tag index lookup in model order
    tag_to_index = {tag_codes_model[i]: i for i in range(len(tag_codes_model))}

    rows = []
    for tag_val in tag_codes_model:
        i_tag = tag_to_index[tag_val]
        t = tag_to_treat.get(tag_val, None)
        if t is None:
            continue
        t = str(t)
        if t not in treat_levels:
            # if df has different string repr, fall back to categorical compare
            t = tag_to_treat[tag_val]

        t_idx = treat_levels.index(t)

        # mu(tag,treat)
        mu_s = mu0 + mu_tag[:, i_tag]
        if mu_treat_full is not None:
            mu_s = mu_s + mu_treat_full[:, t_idx]
        mu_s = circular_modulo(mu_s)

        # kappa(tag,treat)
        eta_k_s = ek0 + (ek_tag[:, i_tag] * sigk)
        if ktreat_full is not None:
            eta_k_s = eta_k_s + ktreat_full[:, t_idx]
        kappa_s = softplus(eta_k_s)

        # p(tag,treat)
        logit_s = lp0 + (lp_tag[:, i_tag] * sigp)
        if beta_full is not None:
            logit_s = logit_s + beta_full[:, t_idx]
        p_s = sigmoid(logit_s)

        # derived
        f0_goal = vonmises_pdf_scalar(0.0, mu_s, kappa_s)
        dens0_s = p_s * f0_goal + (1.0 - p_s) * f0_unif

        A_goal = bessel_ratio(kappa_s)
        align0_s = p_s * (A_goal * np.cos(mu_s)) + (1.0 - p_s) * A_unif

        # summarize helper
        def _summ(x):
            lo, hi = hdi_1d(x, prob=hdi_prob)
            return float(np.mean(x)), lo, hi

        p_m, p_lo, p_hi = _summ(p_s)
        k_m, k_lo, k_hi = _summ(kappa_s)
        d_m, d_lo, d_hi = _summ(dens0_s)
        a_m, a_lo, a_hi = _summ(align0_s)

        rows.append({
            "treatment": t,
            "treatment_label": treat_labels[t_idx],
            "individual": tag_val,

            "p_mean": p_m, "p_lo": p_lo, "p_hi": p_hi,
            "kappa_mean": k_m, "kappa_lo": k_lo, "kappa_hi": k_hi,
            "dens0_mean": d_m, "dens0_lo": d_lo, "dens0_hi": d_hi,
            "align0_mean": a_m, "align0_lo": a_lo, "align0_hi": a_hi,
        })

    out = pd.DataFrame(rows)
    out["treatment"] = pd.Categorical(out["treatment"], categories=treat_levels, ordered=True)
    return out, treat_levels


def pop_summary_by_treatment_from_model_terms(
    idata_mix_full,
    df_fit,
    *,
    treat_col="treatment",
    treat_name_col="treatmentName",
    hdi_prob=0.95,
    max_draws=2000,
    kappa_uniform=KAPPA_UNIFORM,
    seed=123,
):
    post = idata_mix_full.posterior

    for v in ["mu0", "logit_p0", "eta_kappa0"]:
        if v not in post.data_vars:
            raise KeyError(f"Missing '{v}' in posterior. Available: {list(post.data_vars)[:120]} ...")

    treat_levels, labels = _treat_levels_and_labels(df_fit, treat_col=treat_col, treat_name_col=treat_name_col)
    n_treat = len(treat_levels)

    rng = np.random.default_rng(seed)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]).reshape(-1))    # (S,)
    lp0 = _stack_samples_first(post["logit_p0"]).reshape(-1)                # (S,)
    ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)              # (S,)

    def _full_contr(arr_name):
        if arr_name not in post.data_vars:
            return None
        arr = _stack_samples_first(post[arr_name])
        if arr.ndim == 1:
            arr = arr[:, None]
        if arr.shape[1] != (n_treat - 1):
            raise ValueError(f"{arr_name} dim {arr.shape[1]} != n_treat-1={n_treat-1}")
        return np.concatenate([np.zeros((arr.shape[0], 1)), arr], axis=1)   # (S, n_treat)

    mu_treat_full = _full_contr("mu_treat")
    ktreat_full   = _full_contr("kappa_treat")
    beta_full     = _full_contr("beta_treat")

    S = lp0.size
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=max_draws, replace=False)
        mu0 = mu0[idx]; lp0 = lp0[idx]; ek0 = ek0[idx]
        if mu_treat_full is not None: mu_treat_full = mu_treat_full[idx, :]
        if ktreat_full   is not None: ktreat_full   = ktreat_full[idx, :]
        if beta_full     is not None: beta_full     = beta_full[idx, :]
        S = max_draws

    f0_unif = float(vonmises_pdf_scalar(0.0, 0.0, kappa_uniform))
    A_unif  = float(bessel_ratio(kappa_uniform))

    def _summ(x):
        lo, hi = hdi_1d(x, prob=hdi_prob)
        return float(np.mean(x)), lo, hi

    rows = []
    for t_idx, t in enumerate(treat_levels):
        mu_s = mu0.copy()
        if mu_treat_full is not None:
            mu_s = circular_modulo(mu_s + mu_treat_full[:, t_idx])
        else:
            mu_s = circular_modulo(mu_s)

        ek_s = ek0.copy()
        if ktreat_full is not None:
            ek_s = ek_s + ktreat_full[:, t_idx]
        kappa_s = softplus(ek_s)

        lp_s = lp0.copy()
        if beta_full is not None:
            lp_s = lp_s + beta_full[:, t_idx]
        p_s = sigmoid(lp_s)

        f0_goal = vonmises_pdf_scalar(0.0, mu_s, kappa_s)
        dens0_s = p_s * f0_goal + (1.0 - p_s) * f0_unif

        A_goal = bessel_ratio(kappa_s)
        align0_s = p_s * (A_goal * np.cos(mu_s)) + (1.0 - p_s) * A_unif

        p_m, p_lo, p_hi = _summ(p_s)
        k_m, k_lo, k_hi = _summ(kappa_s)
        d_m, d_lo, d_hi = _summ(dens0_s)
        a_m, a_lo, a_hi = _summ(align0_s)

        rows.append({
            "treatment": t,
            "treatment_label": labels[t_idx],

            "p_mean": p_m, "p_lo": p_lo, "p_hi": p_hi,
            "kappa_mean": k_m, "kappa_lo": k_lo, "kappa_hi": k_hi,
            "dens0_mean": d_m, "dens0_lo": d_lo, "dens0_hi": d_hi,
            "align0_mean": a_m, "align0_lo": a_lo, "align0_hi": a_hi,
        })

    out = pd.DataFrame(rows)
    out["treatment"] = pd.Categorical(out["treatment"], categories=treat_levels, ordered=True)
    return out


# ---------------------------------------------------------------------------
# Population-level posterior predictive check (hierarchical bootstrap of
# observed vs simulated per-treatment angle densities)
# ---------------------------------------------------------------------------


@dataclass
class GroupData:
    treat_level: str
    tag_values: np.ndarray
    y_by_tag: List[np.ndarray]



def wrap_angle(x: np.ndarray) -> np.ndarray:
    """Wrap to (-pi, pi]."""
    return (x + np.pi) % (2 * np.pi) - np.pi

def von_mises_kde(data: np.ndarray, points: np.ndarray, kappa: float) -> np.ndarray:
    """
    Von Mises KDE (density per radian). Kernel is VM with concentration kappa.
    """
    data = wrap_angle(np.asarray(data, float))
    points = np.asarray(points, float)
    c = 1.0 / (2.0 * np.pi * np.i0(kappa))
    cos_diffs = np.cos(points[:, None] - data[None, :])
    return c * np.exp(kappa * cos_diffs).mean(axis=1)

def hist_density(data: np.ndarray, bin_edges: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Histogram density (per radian) + bin centers.
    """
    h, _ = np.histogram(data, bins=bin_edges, density=True)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return centers, h

def hdi_1d(x: np.ndarray, prob: float = 0.95) -> Tuple[float, float]:
    """
    1D HDI for a unimodal-ish sample. If you prefer quantiles, swap out.
    """
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    if n == 0:
        return np.nan, np.nan
    x = np.sort(x)
    m = int(np.floor(prob * n))
    if m < 1:
        return np.nan, np.nan
    widths = x[m:] - x[: n - m]
    j = int(np.argmin(widths))
    return float(x[j]), float(x[j + m])

def hdi_along_samples(arr: np.ndarray, prob: float = 0.95) -> Tuple[np.ndarray, np.ndarray]:
    """
    HDI along axis=0 where arr shape = (S, T).
    """
    arr = np.asarray(arr, float)
    S, T = arr.shape
    lo = np.empty(T, float)
    hi = np.empty(T, float)
    for t in range(T):
        lo[t], hi[t] = hdi_1d(arr[:, t], prob=prob)
    return lo, hi


# ============================================================
# Helpers: tqdm + joblib
# ============================================================

@contextmanager
def tqdm_joblib(tqdm_object):
    """
    Context manager to patch joblib to report into tqdm progress bar.
    """
    import joblib
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old
        tqdm_object.close()


# ============================================================
# Grouping utilities
# ============================================================

@dataclass
class GroupData:
    treat_level: str
    tag_values: np.ndarray
    y_by_tag: List[np.ndarray]

def build_groups(
    df: pd.DataFrame,
    y_col: str = "HeadingRadians",
    tag_col: str = "tag",
    treat_col: str = "treatment",
) -> Tuple[List[str], Dict[str, GroupData]]:
    """
    Build per-treatment lists of per-tag y arrays for fast hierarchical bootstrap.
    """
    df = df.copy()

    # ensure treatment categorical (stable ordering)
    treat_cat = pd.Categorical(df[treat_col])
    df[treat_col] = treat_cat
    treat_levels = list(treat_cat.categories)

    groups: Dict[str, GroupData] = {}
    for t in treat_levels:
        sub = df[df[treat_col] == t]
        tags = sub[tag_col].to_numpy()
        y = wrap_angle(sub[y_col].to_numpy())

        uniq_tags = np.unique(tags)
        y_by_tag = []
        for tg in uniq_tags:
            y_by_tag.append(y[tags == tg])
        groups[t] = GroupData(treat_level=str(t), tag_values=uniq_tags, y_by_tag=y_by_tag)

    return treat_levels, groups

def hier_boot_resample_y(y_by_tag: List[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """
    Hierarchical bootstrap:
      - resample tags with replacement
      - within each selected tag, resample observations with replacement (same count as original tag)
    Returns a pooled y* array.
    """
    n_tags = len(y_by_tag)
    tag_idx_star = rng.integers(0, n_tags, size=n_tags)
    out = []
    for j in tag_idx_star:
        yj = y_by_tag[j]
        if yj.size == 0:
            continue
        out.append(rng.choice(yj, size=yj.size, replace=True))
    if not out:
        return np.array([], float)
    return wrap_angle(np.concatenate(out))


# ============================================================
# Observed bootstrap summaries
# ============================================================

def one_obs_boot_rep(
    y_by_tag: List[np.ndarray],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    y_star = hier_boot_resample_y(y_by_tag, rng)
    kde = von_mises_kde(y_star, angles, kappa=kde_kappa)
    _, h = hist_density(y_star, bin_edges)
    return kde, h

def bootstrap_observed(
    groups: Dict[str, GroupData],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    n_boot: int,
    kde_kappa: float,
    n_jobs: int,
    seed: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Returns dict:
      obs_boot[treat]["kde"]  shape (B, P)
      obs_boot[treat]["hist"] shape (B, BINS)
    """
    obs_boot: Dict[str, Dict[str, np.ndarray]] = {}
    base_rng = np.random.default_rng(seed)

    for t, gd in groups.items():
        seeds = base_rng.integers(1, 2**31 - 1, size=n_boot, dtype=np.int64)

        with tqdm_joblib(tqdm(total=n_boot, desc=f"Observed bootstrap ({t})", leave=False)):
            reps = Parallel(n_jobs=n_jobs)(
                delayed(one_obs_boot_rep)(gd.y_by_tag, angles, bin_edges, kde_kappa, int(s))
                for s in seeds
            )

        kde_mat = np.stack([r[0] for r in reps], axis=0)
        hist_mat = np.stack([r[1] for r in reps], axis=0)
        obs_boot[t] = {"kde": kde_mat, "hist": hist_mat}

    return obs_boot


# ============================================================
# Posterior extraction + simulation from fitted model (includes tag REs)
# ============================================================

def _stack_posterior(da) -> np.ndarray:
    """
    Stack chain/draw -> sample and return numpy array.
    """
    da2 = da.stack(sample=("chain", "draw"))
    other = [d for d in da2.dims if d != "sample"]
    da2 = da2.transpose("sample", *other)
    return np.asarray(da2.values)

def extract_simulation_arrays(idata, n_draws: int, seed: int = 123):
    """
    Extract mu, kappa, p from idata.posterior for simulation.

    Expects:
      - mu: shape (S, N)
      - kappa: shape (S, N)
      - p_goal: shape (S, N, 2) or equivalent; we take component 0 as p(directed)
    """
    post = idata.posterior

    if "mu" not in post.data_vars or "kappa" not in post.data_vars or "p_goal" not in post.data_vars:
        raise KeyError(
            "Need posterior vars 'mu', 'kappa', and 'p_goal' in idata. "
            f"Available: {list(post.data_vars)[:80]} ..."
        )

    mu = _stack_posterior(post["mu"])
    kappa = _stack_posterior(post["kappa"])
    p_goal = _stack_posterior(post["p_goal"])

    mu = mu.reshape(mu.shape[0], -1)
    kappa = kappa.reshape(kappa.shape[0], -1)

    if p_goal.ndim == 2:
        p_dir = p_goal.reshape(p_goal.shape[0], -1)
    elif p_goal.ndim == 3:
        if p_goal.shape[-1] == 2:
            p_dir = p_goal[..., 0].reshape(p_goal.shape[0], -1)
        elif p_goal.shape[1] == 2:
            p_dir = p_goal[:, 0, :].reshape(p_goal.shape[0], -1)
        elif p_goal.shape[2] == 2:
            p_dir = p_goal[:, :, 0].reshape(p_goal.shape[0], -1)
        else:
            raise ValueError(f"Unexpected p_goal shape: {p_goal.shape} (need a length-2 axis).")
    else:
        raise ValueError(f"Unexpected p_goal ndim={p_goal.ndim}, shape={p_goal.shape}")

    S = mu.shape[0]
    rng = np.random.default_rng(seed)

    if n_draws is not None and S > n_draws:
        idx = rng.choice(S, size=n_draws, replace=False)
        mu = mu[idx]
        kappa = kappa[idx]
        p_dir = p_dir[idx]

    return mu, kappa, p_dir

def simulate_one_draw(
    mu_s: np.ndarray,
    kappa_s: np.ndarray,
    p_s: np.ndarray,
    rng: np.random.Generator,
    kappa_uniform: float = 1e-6,
) -> np.ndarray:
    """
    Simulate y_rep for ONE posterior draw, for all observations (length N).
    Uses numpy's vonmises sampler.
    """
    N = mu_s.size
    u = rng.random(N)
    is_goal = u < p_s

    y_rep = np.empty(N, float)
    if is_goal.any():
        y_rep[is_goal] = rng.vonmises(mu_s[is_goal], kappa_s[is_goal])
    if (~is_goal).any():
        y_rep[~is_goal] = rng.vonmises(0.0, kappa_uniform, size=(~is_goal).sum())

    return wrap_angle(y_rep)


# ============================================================
# Posterior predictive summaries (simulate -> hierarchical bootstrap summaries)
# ============================================================

def one_ppc_rep_for_treatment(
    y_rep: np.ndarray,
    df: pd.DataFrame,
    treat_level: str,
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    seed: int,
    n_boot_within: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    For a fixed simulated dataset y_rep (length N), compute *n_boot_within* hierarchical-bootstrap
    summaries within one treatment. Returns the average over within-draw boot replicates.
    """
    rng = np.random.default_rng(seed)

    sub = df[df["treatment"] == treat_level]
    idx = sub.index.to_numpy()
    tags = sub["tag"].to_numpy()
    y_t = y_rep[idx]

    uniq_tags = np.unique(tags)
    y_by_tag = [y_t[tags == tg] for tg in uniq_tags]

    kde_acc = []
    hist_acc = []
    for _ in range(n_boot_within):
        y_star = hier_boot_resample_y(y_by_tag, rng)
        kde_acc.append(von_mises_kde(y_star, angles, kappa=kde_kappa))
        _, h = hist_density(y_star, bin_edges)
        hist_acc.append(h)

    return np.mean(np.stack(kde_acc, axis=0), axis=0), np.mean(np.stack(hist_acc, axis=0), axis=0)

def posterior_predictive_summaries(
    idata,
    df: pd.DataFrame,
    treat_levels: List[str],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    n_draws: int,
    n_boot_within_draw: int,
    n_jobs: int,
    seed: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Returns dict:
      ppc[treat]["kde"]  shape (S, P)
      ppc[treat]["hist"] shape (S, BINS)
    where S = n_draws used.
    """
    mu, kappa, p_dir = extract_simulation_arrays(idata, n_draws=n_draws, seed=seed)
    S, N = mu.shape

    if N != len(df):
        raise ValueError(f"Posterior N={N} does not match df length={len(df)}. Check indexing/order.")

    base_rng = np.random.default_rng(seed)
    draw_seeds = base_rng.integers(1, 2**31 - 1, size=S, dtype=np.int64)

    def work_one_draw(s: int):
        rng = np.random.default_rng(int(draw_seeds[s]))
        y_rep = simulate_one_draw(mu[s], kappa[s], p_dir[s], rng, kappa_uniform=KAPPA_UNIFORM)

        out = {}
        treat_seeds = rng.integers(1, 2**31 - 1, size=len(treat_levels), dtype=np.int64)
        for j, t in enumerate(treat_levels):
            kde_j, hist_j = one_ppc_rep_for_treatment(
                y_rep=y_rep,
                df=df,
                treat_level=t,
                angles=angles,
                bin_edges=bin_edges,
                kde_kappa=kde_kappa,
                seed=int(treat_seeds[j]),
                n_boot_within=n_boot_within_draw,
            )
            out[t] = (kde_j, hist_j)
        return out

    with tqdm_joblib(tqdm(total=S, desc="Posterior predictive (draws)", leave=False)):
        draw_results = Parallel(n_jobs=n_jobs)(
            delayed(work_one_draw)(s) for s in range(S)
        )

    ppc: Dict[str, Dict[str, np.ndarray]] = {t: {"kde": [], "hist": []} for t in treat_levels}
    for out in draw_results:
        for t in treat_levels:
            kde_j, hist_j = out[t]
            ppc[t]["kde"].append(kde_j)
            ppc[t]["hist"].append(hist_j)

    for t in treat_levels:
        ppc[t]["kde"] = np.stack(ppc[t]["kde"], axis=0)
        ppc[t]["hist"] = np.stack(ppc[t]["hist"], axis=0)

    return ppc


# ============================================================
# Summaries + plot
# ============================================================

def summarize_band(mat: np.ndarray, prob: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    mat: (S, T)
    returns: median, lo, hi (T,)
    """
    med = np.median(mat, axis=0)
    lo, hi = hdi_along_samples(mat, prob=prob)
    return med, lo, hi

from matplotlib.patches import Circle

def format_polar(ax, edge_radius: float, outer_radius: float):
    ax.spines["polar"].set_visible(False)
    ax.grid(False)
    ax.set_yticks([])
    ax.set_yticklabels([])

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    # Force full circle in [0, 2pi)
    ax.set_xlim(0, 2 * np.pi)

    ax.set_thetagrids(
        [0, 90, 180, 270],
        labels=[" 0°", "", "", ""]
    )

    # keep labels, but remove the tick marks (the little line at the rim)
    ax.tick_params(axis="x", which="major", pad=-2, length=0)

    # Reference circles
    ax.add_artist(Circle((0, 0), edge_radius,
                         transform=ax.transData._b,
                         color="0.35", fill=False, linewidth=0.7))
    ax.add_artist(Circle((0, 0), outer_radius,
                         transform=ax.transData._b,
                         color="0.85", fill=False, linewidth=0.7))

    # Reference ray at 0°
    ax.plot([0, 0], [edge_radius, outer_radius],
            color="0.85", lw=0.8, zorder=0)


def von_mises_kde(data: np.ndarray, points: np.ndarray, kappa: float) -> np.ndarray:
    """
    Von Mises KDE (density per radian). Kernel is VM with concentration kappa.
    """
    data = wrap_angle(np.asarray(data, float))
    points = np.asarray(points, float)
    c = 1.0 / (2.0 * np.pi * np.i0(kappa))
    cos_diffs = np.cos(points[:, None] - data[None, :])
    return c * np.exp(kappa * cos_diffs).mean(axis=1)


def hist_density(data: np.ndarray, bin_edges: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Histogram density (per radian) + bin centers.
    """
    h, _ = np.histogram(data, bins=bin_edges, density=True)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return centers, h


def tqdm_joblib(tqdm_object):
    """
    Context manager to patch joblib to report into tqdm progress bar.
    """
    import joblib
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old
        tqdm_object.close()


def build_groups(
    df: pd.DataFrame,
    y_col: str = "HeadingRadians",
    tag_col: str = "tag",
    treat_col: str = "treatment",
) -> Tuple[List[str], Dict[str, GroupData]]:
    """
    Build per-treatment lists of per-tag y arrays for fast hierarchical bootstrap.
    """
    df = df.copy()

    # ensure treatment categorical (stable ordering)
    treat_cat = pd.Categorical(df[treat_col])
    df[treat_col] = treat_cat
    treat_levels = list(treat_cat.categories)

    groups: Dict[str, GroupData] = {}
    for t in treat_levels:
        sub = df[df[treat_col] == t]
        tags = sub[tag_col].to_numpy()
        y = wrap_angle(sub[y_col].to_numpy())

        uniq_tags = np.unique(tags)
        y_by_tag = []
        for tg in uniq_tags:
            y_by_tag.append(y[tags == tg])
        groups[t] = GroupData(treat_level=str(t), tag_values=uniq_tags, y_by_tag=y_by_tag)

    return treat_levels, groups


def hier_boot_resample_y(y_by_tag: List[np.ndarray], rng: np.random.Generator) -> np.ndarray:
    """
    Hierarchical bootstrap:
      - resample tags with replacement
      - within each selected tag, resample observations with replacement (same count as original tag)
    Returns a pooled y* array.
    """
    n_tags = len(y_by_tag)
    tag_idx_star = rng.integers(0, n_tags, size=n_tags)
    out = []
    for j in tag_idx_star:
        yj = y_by_tag[j]
        if yj.size == 0:
            continue
        out.append(rng.choice(yj, size=yj.size, replace=True))
    if not out:
        return np.array([], float)
    return wrap_angle(np.concatenate(out))


def one_obs_boot_rep(
    y_by_tag: List[np.ndarray],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    y_star = hier_boot_resample_y(y_by_tag, rng)
    kde = von_mises_kde(y_star, angles, kappa=kde_kappa)
    _, h = hist_density(y_star, bin_edges)
    return kde, h


def bootstrap_observed(
    groups: Dict[str, GroupData],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    n_boot: int,
    kde_kappa: float,
    n_jobs: int,
    seed: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Returns dict:
      obs_boot[treat]["kde"]  shape (B, P)
      obs_boot[treat]["hist"] shape (B, BINS)
    """
    obs_boot: Dict[str, Dict[str, np.ndarray]] = {}
    base_rng = np.random.default_rng(seed)

    for t, gd in groups.items():
        seeds = base_rng.integers(1, 2**31 - 1, size=n_boot, dtype=np.int64)

        with tqdm_joblib(tqdm(total=n_boot, desc=f"Observed bootstrap ({t})", leave=False)):
            reps = Parallel(n_jobs=n_jobs)(
                delayed(one_obs_boot_rep)(gd.y_by_tag, angles, bin_edges, kde_kappa, int(s))
                for s in seeds
            )

        kde_mat = np.stack([r[0] for r in reps], axis=0)
        hist_mat = np.stack([r[1] for r in reps], axis=0)
        obs_boot[t] = {"kde": kde_mat, "hist": hist_mat}

    return obs_boot


def _stack_posterior(da) -> np.ndarray:
    """
    Stack chain/draw -> sample and return numpy array.
    """
    da2 = da.stack(sample=("chain", "draw"))
    other = [d for d in da2.dims if d != "sample"]
    da2 = da2.transpose("sample", *other)
    return np.asarray(da2.values)


def extract_simulation_arrays(idata, n_draws: int, seed: int = 123):
    """
    Extract mu, kappa, p from idata.posterior for simulation.

    Expects:
      - mu: shape (S, N)
      - kappa: shape (S, N)
      - p_goal: shape (S, N, 2) or equivalent; we take component 0 as p(directed)
    """
    post = idata.posterior

    if "mu" not in post.data_vars or "kappa" not in post.data_vars or "p_goal" not in post.data_vars:
        raise KeyError(
            "Need posterior vars 'mu', 'kappa', and 'p_goal' in idata. "
            f"Available: {list(post.data_vars)[:80]} ..."
        )

    mu = _stack_posterior(post["mu"])
    kappa = _stack_posterior(post["kappa"])
    p_goal = _stack_posterior(post["p_goal"])

    mu = mu.reshape(mu.shape[0], -1)
    kappa = kappa.reshape(kappa.shape[0], -1)

    if p_goal.ndim == 2:
        p_dir = p_goal.reshape(p_goal.shape[0], -1)
    elif p_goal.ndim == 3:
        if p_goal.shape[-1] == 2:
            p_dir = p_goal[..., 0].reshape(p_goal.shape[0], -1)
        elif p_goal.shape[1] == 2:
            p_dir = p_goal[:, 0, :].reshape(p_goal.shape[0], -1)
        elif p_goal.shape[2] == 2:
            p_dir = p_goal[:, :, 0].reshape(p_goal.shape[0], -1)
        else:
            raise ValueError(f"Unexpected p_goal shape: {p_goal.shape} (need a length-2 axis).")
    else:
        raise ValueError(f"Unexpected p_goal ndim={p_goal.ndim}, shape={p_goal.shape}")

    S = mu.shape[0]
    rng = np.random.default_rng(seed)

    if n_draws is not None and S > n_draws:
        idx = rng.choice(S, size=n_draws, replace=False)
        mu = mu[idx]
        kappa = kappa[idx]
        p_dir = p_dir[idx]

    return mu, kappa, p_dir


def simulate_one_draw(
    mu_s: np.ndarray,
    kappa_s: np.ndarray,
    p_s: np.ndarray,
    rng: np.random.Generator,
    kappa_uniform: float = 1e-6,
) -> np.ndarray:
    """
    Simulate y_rep for ONE posterior draw, for all observations (length N).
    Uses numpy's vonmises sampler.
    """
    N = mu_s.size
    u = rng.random(N)
    is_goal = u < p_s

    y_rep = np.empty(N, float)
    if is_goal.any():
        y_rep[is_goal] = rng.vonmises(mu_s[is_goal], kappa_s[is_goal])
    if (~is_goal).any():
        y_rep[~is_goal] = rng.vonmises(0.0, kappa_uniform, size=(~is_goal).sum())

    return wrap_angle(y_rep)


def one_ppc_rep_for_treatment(
    y_rep: np.ndarray,
    df: pd.DataFrame,
    treat_level: str,
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    seed: int,
    n_boot_within: int = 1,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    For a fixed simulated dataset y_rep (length N), compute *n_boot_within* hierarchical-bootstrap
    summaries within one treatment. Returns the average over within-draw boot replicates.
    """
    rng = np.random.default_rng(seed)

    sub = df[df["treatment"] == treat_level]
    idx = sub.index.to_numpy()
    tags = sub["tag"].to_numpy()
    y_t = y_rep[idx]

    uniq_tags = np.unique(tags)
    y_by_tag = [y_t[tags == tg] for tg in uniq_tags]

    kde_acc = []
    hist_acc = []
    for _ in range(n_boot_within):
        y_star = hier_boot_resample_y(y_by_tag, rng)
        kde_acc.append(von_mises_kde(y_star, angles, kappa=kde_kappa))
        _, h = hist_density(y_star, bin_edges)
        hist_acc.append(h)

    return np.mean(np.stack(kde_acc, axis=0), axis=0), np.mean(np.stack(hist_acc, axis=0), axis=0)


def posterior_predictive_summaries(
    idata,
    df: pd.DataFrame,
    treat_levels: List[str],
    angles: np.ndarray,
    bin_edges: np.ndarray,
    kde_kappa: float,
    n_draws: int,
    n_boot_within_draw: int,
    n_jobs: int,
    seed: int,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Returns dict:
      ppc[treat]["kde"]  shape (S, P)
      ppc[treat]["hist"] shape (S, BINS)
    where S = n_draws used.
    """
    mu, kappa, p_dir = extract_simulation_arrays(idata, n_draws=n_draws, seed=seed)
    S, N = mu.shape

    if N != len(df):
        raise ValueError(f"Posterior N={N} does not match df length={len(df)}. Check indexing/order.")

    base_rng = np.random.default_rng(seed)
    draw_seeds = base_rng.integers(1, 2**31 - 1, size=S, dtype=np.int64)

    def work_one_draw(s: int):
        rng = np.random.default_rng(int(draw_seeds[s]))
        y_rep = simulate_one_draw(mu[s], kappa[s], p_dir[s], rng, kappa_uniform=KAPPA_UNIFORM)

        out = {}
        treat_seeds = rng.integers(1, 2**31 - 1, size=len(treat_levels), dtype=np.int64)
        for j, t in enumerate(treat_levels):
            kde_j, hist_j = one_ppc_rep_for_treatment(
                y_rep=y_rep,
                df=df,
                treat_level=t,
                angles=angles,
                bin_edges=bin_edges,
                kde_kappa=kde_kappa,
                seed=int(treat_seeds[j]),
                n_boot_within=n_boot_within_draw,
            )
            out[t] = (kde_j, hist_j)
        return out

    with tqdm_joblib(tqdm(total=S, desc="Posterior predictive (draws)", leave=False)):
        draw_results = Parallel(n_jobs=n_jobs)(
            delayed(work_one_draw)(s) for s in range(S)
        )

    ppc: Dict[str, Dict[str, np.ndarray]] = {t: {"kde": [], "hist": []} for t in treat_levels}
    for out in draw_results:
        for t in treat_levels:
            kde_j, hist_j = out[t]
            ppc[t]["kde"].append(kde_j)
            ppc[t]["hist"].append(hist_j)

    for t in treat_levels:
        ppc[t]["kde"] = np.stack(ppc[t]["kde"], axis=0)
        ppc[t]["hist"] = np.stack(ppc[t]["hist"], axis=0)

    return ppc


def format_treatment_title(code, label_map, drop_control_if_mixed=True):
    raw = str(label_map.get(code, code)).strip()
    if not raw:
        return str(code)

    parts = [p.strip() for p in raw.split(",") if p.strip()]

    if drop_control_if_mixed and len(parts) > 1:
        parts = [p for p in parts if p.lower() != "control"]

    parts = [p.lower().title() for p in parts]

    return ", ".join(parts) if parts else str(code)


# ---------------------------------------------------------------------------
# Individual-level posterior predictive check (per-tag bootstrap KDE vs
# simulated KDE, one panel per tag)
# ---------------------------------------------------------------------------

RNG_SEED = 123
EDGE_RADIUS = 0.55    # radial offset for the polar panels
AMP = 1.00            # polar radial gain
THETA_N = 360         # angles per panel
HDI_PROB = 0.95
KAPPA_KDE = 15.0      # von Mises KDE concentration for the per-tag bands
N_BOOT = 500          # observed-bootstrap resamples per tag
MAX_DRAWS = 500       # posterior subsample size per tag
N_JOBS = -1
PREFER = "threads"


def individual_ppc(idata_mix_full, df_track, *, n_cols=8):
    """Per-tag PPC: observed hierarchical-style bootstrap KDE vs the
    posterior-predictive KDE of each tag's goal-relative headings, with
    per-tag median + 95% HDI bands rescaled to the polar radius.

    Subsamples the vMU posterior to ``MAX_DRAWS`` draws and evaluates
    every tag in parallel (threads; the von Mises draws run in numpy).
    Returns a dict consumed by
    ``goldshtein_bats_plots.plot_individual_ppc_grid``: ``results``
    (one row per tag with the band arrays), ``angles``, ``to_r``
    (density -> polar radius), ``RMAX``, the treatment labels/colors and
    the tag ordering.
    """
    df = df_track.copy()

    treat_cat = pd.Categorical(df["treatment"])
    df["treatment"] = treat_cat
    treat_levels = list(treat_cat.categories)
    n_treat = len(treat_levels)

    if "treatmentName" in df.columns:
        treat_name_map = df.groupby("treatment", observed=False)["treatmentName"].first().to_dict()
    else:
        treat_name_map = {t: str(t) for t in treat_levels}

    nice_treat_label = {t: format_treatment_nicely(treat_name_map.get(t, t)) for t in treat_levels}

    tag_col   = "tag"
    y_col     = "HeadingRadians"
    treat_col = "treatment"

    tag_to_rows  = df.groupby(tag_col, observed=False).groups
    tag_to_treat = df.groupby(tag_col, observed=False)[treat_col].first().to_dict()

    tags = list(tag_to_rows.keys())

    # Colors by treatment (used later)
    cmap = mpl.colormaps.get_cmap("tab10")
    treat_to_color = {t: cmap(i % 10) for i, t in enumerate(treat_levels)}

    # Evaluation grid: full circle [0, 2π)
    angles = np.linspace(0.0, 2.0 * np.pi, THETA_N, endpoint=False)

    # ============================================================
    # Extract posterior arrays needed to compute per-tag parameters
    # ============================================================
    post = idata_mix_full.posterior

    required = ["mu0", "mu_tag", "eta_kappa0", "eta_kappa_tag", "sigma_kappa_tag",
                "logit_p0", "logit_p_tag", "sigma_p_tag"]
    missing = [v for v in required if v not in post.data_vars]
    if missing:
        raise KeyError(f"Missing {missing}. Available: {list(post.data_vars)[:80]} ...")

    has_mu_treat    = ("mu_treat" in post.data_vars) and (n_treat > 1)
    has_kappa_treat = ("kappa_treat" in post.data_vars) and (n_treat > 1)
    has_beta_treat  = ("beta_treat" in post.data_vars) and (n_treat > 1)

    mu0 = circular_modulo(_stack_samples_first(post["mu0"]).reshape(-1))
    mu_tag = _stack_samples_first(post["mu_tag"])
    eta_kappa0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)
    eta_kappa_tag = _stack_samples_first(post["eta_kappa_tag"])
    sigma_kappa_tag = _stack_samples_first(post["sigma_kappa_tag"]).reshape(-1)

    logit_p0 = _stack_samples_first(post["logit_p0"]).reshape(-1)
    logit_p_tag = _stack_samples_first(post["logit_p_tag"])
    sigma_p_tag = _stack_samples_first(post["sigma_p_tag"]).reshape(-1)

    if has_mu_treat:
        mu_treat = _stack_samples_first(post["mu_treat"])
        if mu_treat.ndim == 1:
            mu_treat = mu_treat[:, None]
    if has_kappa_treat:
        kappa_treat = _stack_samples_first(post["kappa_treat"])
        if kappa_treat.ndim == 1:
            kappa_treat = kappa_treat[:, None]
    if has_beta_treat:
        beta_treat = _stack_samples_first(post["beta_treat"])
        if beta_treat.ndim == 1:
            beta_treat = beta_treat[:, None]

    # Subsample posterior draws for PPC
    rng_global = np.random.default_rng(RNG_SEED)
    S = mu0.shape[0]
    if (MAX_DRAWS is not None) and (S > MAX_DRAWS):
        idx = rng_global.choice(S, size=MAX_DRAWS, replace=False)
        mu0 = mu0[idx]
        mu_tag = mu_tag[idx, :]
        eta_kappa0 = eta_kappa0[idx]
        eta_kappa_tag = eta_kappa_tag[idx, :]
        sigma_kappa_tag = sigma_kappa_tag[idx]
        logit_p0 = logit_p0[idx]
        logit_p_tag = logit_p_tag[idx, :]
        sigma_p_tag = sigma_p_tag[idx]
        if has_mu_treat:
            mu_treat = mu_treat[idx, :]
        if has_kappa_treat:
            kappa_treat = kappa_treat[idx, :]
        if has_beta_treat:
            beta_treat = beta_treat[idx, :]
        S = MAX_DRAWS

    tag_codes_model = np.unique(df[tag_col].to_numpy())
    tag_to_model_idx = {tag: i for i, tag in enumerate(tag_codes_model)}

    # ============================================================
    # Compute one individual's PPC (returns KDE summaries)
    # ============================================================
    def compute_one_tag(tag_value, seed=0):
        rng = np.random.default_rng(seed)

        rows = tag_to_rows[tag_value]
        y_obs = df.loc[rows, y_col].to_numpy(dtype=float)
        y_obs = wrap_0_2pi(y_obs)

        n_obs = y_obs.size
        treat = tag_to_treat[tag_value]

        if n_obs < 5:
            obs_mat = np.full((1, THETA_N), np.nan)
            ppc_mat = np.full((1, THETA_N), np.nan)
            obs_med, obs_lo, obs_hi = quantile_band(obs_mat, prob=HDI_PROB)
            ppc_med, ppc_lo, ppc_hi = quantile_band(ppc_mat, prob=HDI_PROB)
            return dict(
                tag=tag_value, treat=treat,
                obs_med=obs_med, obs_lo=obs_lo, obs_hi=obs_hi,
                ppc_med=ppc_med, ppc_lo=ppc_lo, ppc_hi=ppc_hi,
                peak=float(np.nanmax([np.nanmax(obs_hi), np.nanmax(ppc_hi)]))
            )

        treat_code = treat_levels.index(treat)

        if tag_value not in tag_to_model_idx:
            raise KeyError(f"Tag {tag_value} not found in model tag mapping.")
        j_tag = tag_to_model_idx[tag_value]

        mu_s = mu0 + mu_tag[:, j_tag]
        if has_mu_treat and (treat_code > 0):
            mu_s = mu_s + mu_treat[:, treat_code - 1]
        mu_s = circular_modulo(mu_s)

        eta_k_s = eta_kappa0 + eta_kappa_tag[:, j_tag] * sigma_kappa_tag
        if has_kappa_treat and (treat_code > 0):
            eta_k_s = eta_k_s + kappa_treat[:, treat_code - 1]
        kappa_s = softplus(eta_k_s)

        logit_p_s = logit_p0 + logit_p_tag[:, j_tag] * sigma_p_tag
        if has_beta_treat and (treat_code > 0):
            logit_p_s = logit_p_s + beta_treat[:, treat_code - 1]
        p_goal_s = sigmoid(logit_p_s)

        # Observed bootstrap KDE
        obs_kdes = np.empty((N_BOOT, THETA_N), dtype=float)
        for b in range(N_BOOT):
            boot = y_obs[rng.integers(0, n_obs, size=n_obs)]
            obs_kdes[b, :] = von_mises_kde(boot, angles, kappa=KAPPA_KDE)
        obs_med, obs_lo, obs_hi = quantile_band(obs_kdes, prob=HDI_PROB)

        # Posterior predictive KDE
        ppc_kdes = np.empty((S, THETA_N), dtype=float)
        for s in range(S):
            u = rng.uniform(size=n_obs)
            is_goal = (u < p_goal_s[s])

            n_goal = int(is_goal.sum())
            n_unif = n_obs - n_goal

            sim = np.empty(n_obs, dtype=float)
            if n_goal > 0:
                sim[is_goal] = rng.vonmises(mu_s[s], kappa_s[s], size=n_goal)
            if n_unif > 0:
                sim[~is_goal] = rng.vonmises(0.0, KAPPA_UNIFORM, size=n_unif)

            sim = wrap_0_2pi(sim)
            ppc_kdes[s, :] = von_mises_kde(sim, angles, kappa=KAPPA_KDE)

        ppc_med, ppc_lo, ppc_hi = quantile_band(ppc_kdes, prob=HDI_PROB)
        peak = float(np.nanmax([obs_hi.max(), ppc_hi.max()]))

        return dict(
            tag=tag_value, treat=treat,
            obs_med=obs_med, obs_lo=obs_lo, obs_hi=obs_hi,
            ppc_med=ppc_med, ppc_lo=ppc_lo, ppc_hi=ppc_hi,
            peak=peak
        )

    # ============================================================
    # Run compute (results used by plotting cell)
    # ============================================================
    results = Parallel(n_jobs=N_JOBS, prefer=PREFER)(
        delayed(compute_one_tag)(tag, seed=10_000 + i)
        for i, tag in enumerate(tqdm(tags, desc="Individual PPC (bootstrap KDE + PPC KDE)"))
    )

    treat_order = {t: i for i, t in enumerate(treat_levels)}
    results = sorted(results, key=lambda r: treat_order.get(r["treat"], 10**9))

    peak_max = max(r["peak"] for r in results if np.isfinite(r["peak"])) + 1e-12

    def to_r(density):
        return EDGE_RADIUS + AMP * (density / peak_max)

    RMAX = EDGE_RADIUS + AMP * 1.01

    return {
        "results": results,
        "angles": angles,
        "to_r": to_r,
        "RMAX": RMAX,
        "EDGE_RADIUS": EDGE_RADIUS,
        "AMP": AMP,
        "treat_levels": treat_levels,
        "treat_to_color": treat_to_color,
        "nice_treat_label": nice_treat_label,
        "treat_order": treat_order,
    }

# ---------------------------------------------------------------------------
# Per-step posterior responsibility p(goal component | step)
# ---------------------------------------------------------------------------


# The first (logit-space) implementation, summarized on the logit scale
# (the "summarize_on" quirk of the annotated notebook is kept verbatim);
# it drives the p_goal maps of the per-treatment track panels.

def posterior_prob_goal_per_step(
    idata,
    y_obs,
    *,
    max_draws=1200,
    hdi_prob=0.95,
    eps=1e-9,
    summarize_on="logit",  # "logit" or "prob"
):
    """
    Posterior responsibility for a 2-component mixture:
      goal ~ von Mises(mu, kappa)
      other ~ UniformCircle (density 1/(2π))

    r = p(goal | y, draw) = sigmoid( logit(w) + log f_goal(y) - log f_unif(y) )

    With proper uniform on the circle:
      log f_goal - log f_unif = kappa*cos(y-mu) - log I0(kappa)
    """
    post = idata.posterior
    y_obs = np.asarray(y_obs, float).reshape(-1)
    n_obs = y_obs.size

    # --- mixture weight w(goal) per obs per draw ---
    if "p_goal" in post.data_vars:
        w_full = _as_sample_obs(post["p_goal"], n_obs, want_components=True)  # (S, n_obs, 2)
        w = w_full[:, :, 0]
    else:
        if "logit_p" not in post.data_vars:
            raise KeyError("Need either posterior var 'p_goal' (preferred) or 'logit_p'.")
        logit_p = _as_sample_obs(post["logit_p"], n_obs, want_components=False)
        w = sps.expit(logit_p)

    # --- mu, kappa per obs per draw ---
    if "mu" in post.data_vars:
        mu = _as_sample_obs(post["mu"], n_obs, want_components=False)
    elif "mu0" in post.data_vars:
        mu0 = _stack_samples_first(post["mu0"]).values.astype(float).reshape(-1)  # (S,)
        mu = np.broadcast_to(mu0[:, None], (mu0.size, n_obs))
    else:
        raise KeyError("Need posterior var 'mu' or 'mu0'.")

    if "kappa" in post.data_vars:
        kappa = _as_sample_obs(post["kappa"], n_obs, want_components=False)
    elif "eta_kappa0" in post.data_vars:
        ek0 = _stack_samples_first(post["eta_kappa0"]).values.astype(float).reshape(-1)  # (S,)
        k0 = _softplus(ek0)
        kappa = np.broadcast_to(k0[:, None], (k0.size, n_obs))
    else:
        raise KeyError("Need posterior var 'kappa' or 'eta_kappa0'.")

    # --- subsample draws once ---
    S = w.shape[0]
    rng = np.random.default_rng(123)
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=max_draws, replace=False)
        w = w[idx, :]
        mu = mu[idx, :]
        kappa = kappa[idx, :]
        S = max_draws

    # --- compute responsibility logit using exact uniform circle ---
    w = np.clip(w, eps, 1.0 - eps)
    logit_w = np.log(w) - np.log1p(-w)  # (S, n_obs)

    y = y_obs[None, :]
    ll_diff = kappa * np.cos(y - mu) - _log_i0(kappa)  # log f_VM - log f_unif

    logit_r = logit_w + ll_diff
    resp = sps.expit(logit_r)

    if summarize_on == "prob":
        mean_r = resp.mean(axis=0)
        hdi = az.hdi(resp, prob=hdi_prob, skipna=True)
        return mean_r, hdi, resp

    if summarize_on == "logit":
        mean_logit = logit_r.mean(axis=0)
        hdi_logit = az.hdi(logit_r, prob=hdi_prob, skipna=True)

        mean_r = sps.expit(mean_logit * 2)
        hdi = sps.expit(hdi_logit)
        return mean_r, hdi, resp

    raise ValueError("summarize_on must be 'prob' or 'logit'")


# The curves implementation: draws the per-step responsibility and
# smooths it along the trajectory into population / per-tag kernel
# curves over distance-to-destination.

def _next_pow_two(n: int) -> int:
    return 1 << (int(n) - 1).bit_length()


def _fft_convolve_same_1d(x, k):
    """1D convolution (same) using FFT. x,k are 1D arrays."""
    x = np.asarray(x, float)
    k = np.asarray(k, float)
    nx = x.size
    nk = k.size
    n = nx + nk - 1
    nfft = _next_pow_two(n)

    X = np.fft.rfft(x, nfft)
    K = np.fft.rfft(k, nfft)
    y = np.fft.irfft(X * K, nfft)[:n]

    start = (nk - 1) // 2
    return y[start:start + nx]


def _fft_convolve_same_rows(A, k):
    """
    Row-wise convolution (same) using FFT.
    A: (S,G), k: (K,)
    returns: (S,G)
    """
    A = np.asarray(A, float)
    k = np.asarray(k, float)
    S, G = A.shape
    K = k.size
    n = G + K - 1
    nfft = _next_pow_two(n)

    Af = np.fft.rfft(A, nfft, axis=1)
    kf = np.fft.rfft(k, nfft)
    y = np.fft.irfft(Af * kf[None, :], nfft, axis=1)[:, :n]

    start = (K - 1) // 2
    return y[:, start:start + G]


def weighted_gaussian_smooth_draws(draw_curve_SxG, xg_G, weight_G, bw_x, eps=1e-12):
    """
    Smooth each draw curve along x with weighted gaussian smoothing:
      y_sm(x) = conv(y*w, k) / conv(w, k)
    where w is (weight_G) and k is gaussian kernel on grid.
    Does not plot beyond support if you mask afterward.
    """
    y = np.asarray(draw_curve_SxG, float)  # (S,G)
    xg = np.asarray(xg_G, float)
    w = np.asarray(weight_G, float)

    if y.ndim != 2 or xg.ndim != 1 or w.ndim != 1:
        raise ValueError("Bad shapes for weighted_gaussian_smooth_draws.")
    if y.shape[1] != xg.size or w.size != xg.size:
        raise ValueError("draw_curve, xg, weight must share grid length.")

    if (bw_x is None) or (not np.isfinite(bw_x)) or (bw_x <= 0):
        return y.copy()

    # uniform grid spacing assumed
    dx = float(np.nanmedian(np.diff(xg)))
    if not np.isfinite(dx) or dx <= 0:
        return y.copy()

    sigma = float(bw_x) / dx
    if sigma < 0.5:  # extremely narrow -> not helpful
        return y.copy()

    radius = int(np.ceil(4.0 * sigma))
    idx = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-0.5 * (idx / sigma) ** 2)

    # numerator: y*w (NaNs treated as 0 with weight 0)
    y0 = np.where(np.isfinite(y), y, 0.0)
    num = y0 * w[None, :]

    # convolve numerator rows + denominator once
    num_sm = _fft_convolve_same_rows(num, k)
    den_sm = _fft_convolve_same_1d(w, k)

    out = num_sm / (den_sm[None, :] + eps)
    out[:, den_sm < eps] = np.nan
    return out


# the per-curve version of the responsibility (returns the draws)

def _posterior_prob_goal_per_step_curves(
    idata,
    y_obs,
    *,
    kappa_uniform=1e-6,
    goal_comp_index=0,
    max_draws=1200,
    seed=123,
    hdi_prob=0.95,
):
    post = idata.posterior
    y_obs = np.asarray(y_obs, float)
    n_obs = y_obs.size
    rng = np.random.default_rng(seed)

    if "p_goal" not in post.data_vars:
        raise KeyError("Expected deterministic 'p_goal' in posterior.")

    w = _stack_samples_first(post["p_goal"])  # (S, ...)

    if w.ndim == 2:
        axes = list(range(1, w.ndim))
        obs_axes = [ax for ax in axes if w.shape[ax] == n_obs]
        if len(obs_axes) != 1:
            if w.shape[1] != n_obs:
                raise ValueError(f"Can't identify obs axis for p_goal. shape={w.shape}, n_obs={n_obs}")
            w_goal = w
        else:
            w_goal = np.moveaxis(w, obs_axes[0], 1).reshape(w.shape[0], n_obs)

    elif w.ndim >= 3:
        axes = list(range(1, w.ndim))
        obs_axes = [ax for ax in axes if w.shape[ax] == n_obs]
        comp_axes = [ax for ax in axes if w.shape[ax] == 2]
        if len(obs_axes) != 1 or len(comp_axes) != 1:
            raise ValueError(f"Can't identify obs/comp dims for p_goal. shape={w.shape}, n_obs={n_obs}")

        w = np.moveaxis(w, obs_axes[0], 1)
        comp_candidates = [ax for ax in range(2, w.ndim) if w.shape[ax] == 2]
        if len(comp_candidates) != 1:
            raise ValueError(f"Can't identify comp axis after moving obs. shape={w.shape}")
        w = np.moveaxis(w, comp_candidates[0], 2)  # (S, n_obs, 2)

        if goal_comp_index not in (0, 1):
            raise ValueError("goal_comp_index must be 0 or 1 when comp dim exists.")
        w_goal = w[:, :, int(goal_comp_index)]
    else:
        raise ValueError(f"Unexpected p_goal shape={w.shape}")

    w_goal = np.clip(w_goal, 1e-12, 1 - 1e-12)  # (S, n_obs)

    # mu
    if "mu" in post.data_vars:
        mu = _stack_samples_first(post["mu"])
        axes = list(range(1, mu.ndim))
        obs_axes = [ax for ax in axes if mu.shape[ax] == n_obs]
        if len(obs_axes) != 1:
            raise ValueError(f"Can't find obs axis in 'mu'. shape={mu.shape}, n_obs={n_obs}")
        mu = np.moveaxis(mu, obs_axes[0], 1).reshape(mu.shape[0], n_obs)
        mu = circular_modulo(mu)
    elif "mu0" in post.data_vars:
        mu0 = _stack_samples_first(post["mu0"]).reshape(-1)
        mu = circular_modulo(mu0[:, None] + np.zeros((mu0.size, n_obs)))
    else:
        raise KeyError("Need either 'mu' or 'mu0' in posterior.")

    # kappa
    if "kappa" in post.data_vars:
        kappa = _stack_samples_first(post["kappa"])
        axes = list(range(1, kappa.ndim))
        obs_axes = [ax for ax in axes if kappa.shape[ax] == n_obs]
        if len(obs_axes) != 1:
            raise ValueError(f"Can't find obs axis in 'kappa'. shape={kappa.shape}, n_obs={n_obs}")
        kappa = np.moveaxis(kappa, obs_axes[0], 1).reshape(kappa.shape[0], n_obs)
        kappa = np.clip(kappa, 1e-12, np.inf)
    elif "eta_kappa0" in post.data_vars:
        ek0 = _stack_samples_first(post["eta_kappa0"]).reshape(-1)
        k0 = softplus(ek0)
        kappa = k0[:, None] + np.zeros((k0.size, n_obs))
    else:
        raise KeyError("Need either 'kappa' or 'eta_kappa0' in posterior.")

    # subsample draws
    S = w_goal.shape[0]
    if (max_draws is not None) and (S > max_draws):
        idx = rng.choice(S, size=int(max_draws), replace=False)
        w_goal = w_goal[idx]
        mu = mu[idx]
        kappa = kappa[idx]

    logit_w = np.log(w_goal) - np.log1p(-w_goal)
    log_f_goal = log_vonmises_pdf(y_obs[None, :], mu, kappa)
    log_f_unif = log_vonmises_pdf(y_obs[None, :], 0.0, kappa_uniform)

    logit_r = logit_w + (log_f_goal - log_f_unif)
    resp_draws = sigmoid(logit_r)

    resp_mean = resp_draws.mean(axis=0)
    resp_hdi = hdi_along_samples_2col(resp_draws, prob=hdi_prob)

    return resp_draws, resp_mean, resp_hdi


def append_pgoal_post_to_df(
    df_track,
    idata,
    *,
    y_col="HeadingRadians",
    kappa_uniform=1e-6,
    goal_comp_index=0,
    max_draws=1200,
    hdi_prob=0.95,
):
    y_obs = pd.to_numeric(df_track[y_col], errors="coerce").to_numpy(float)
    resp_draws, resp_mean, resp_hdi = _posterior_prob_goal_per_step_curves(
        idata,
        y_obs,
        kappa_uniform=kappa_uniform,
        goal_comp_index=goal_comp_index,
        max_draws=max_draws,
        hdi_prob=hdi_prob
    )
    out = df_track.copy()
    out["p_goal_post_mean"] = resp_mean
    out["p_goal_post_lo"] = resp_hdi[:, 0]
    out["p_goal_post_hi"] = resp_hdi[:, 1]
    return out, resp_draws


def population_curve_kernel_per_tag(
    df,
    resp_draws,
    *,
    tag_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",

    x_kind="distance",                 # "distance" or "steps_from_end"
    last_n_steps=150,
    steps_x_units="m",                 # "m" or "steps"
    step_length_m=60.0,

    dist_col="dist_to_dest_m",
    dist_clip_q=None,

    n_grid=2000,
    bandwidth=240.0,
    chunk_size=4000,

    hdi_prob=0.95,
    eps=1e-12,

    # truncation controls
    truncate_to_data_range=True,
    truncate_by_weight=True,
    den_frac=0.01,

    # population stabilization
    pop_post_smooth=True,
    pop_smooth_bandwidth=None,
    min_tags_at_x=1,

    return_tag_summaries=False,
):
    d = df.copy().reset_index(drop=True)
    n_obs = len(d)

    resp = np.asarray(resp_draws, float)
    if resp.ndim != 2 or resp.shape[1] != n_obs:
        raise ValueError(f"resp_draws must be (S, n_obs) matching df rows. Got {resp.shape} vs n_obs={n_obs}")

    d["_rowpos"] = np.arange(n_obs, dtype=int)

    # treatment order + labels
    if isinstance(d[treat_col].dtype, pd.CategoricalDtype):
        treat_levels = list(d[treat_col].cat.categories)
    else:
        treat_levels = sorted(d[treat_col].dropna().unique().tolist())

    if treat_name_col in d.columns:
        name_map = d.groupby(treat_col, observed=False)[treat_name_col].first().to_dict()
        labels = {t: _format_treatment_name(name_map.get(t, t)) for t in treat_levels}
    else:
        labels = {t: _format_treatment_name(t) for t in treat_levels}

    # ---- define x ----
    if x_kind == "steps_from_end":
        d["_step_fwd"] = d.groupby(tag_col, observed=False).cumcount()
        d["_last_fwd"] = d.groupby(tag_col, observed=False)["_step_fwd"].transform("max")
        d["_x_steps"] = d["_step_fwd"] - d["_last_fwd"]  # last=0; earlier negative

        d = d.loc[(d["_x_steps"] >= -int(last_n_steps)) & (d["_x_steps"] <= 0)].copy()

        if str(steps_x_units).lower() in ("m", "meter", "meters"):
            d["_x"] = d["_x_steps"].astype(float) * float(step_length_m)
        elif str(steps_x_units).lower() in ("step", "steps"):
            d["_x"] = d["_x_steps"].astype(float)
        else:
            raise ValueError("steps_x_units must be 'm' or 'steps'")

    elif x_kind == "distance":
        if dist_col not in d.columns:
            raise KeyError(f"Missing distance column '{dist_col}' in df.")
        xvals = pd.to_numeric(d[dist_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(xvals)
        d = d.loc[ok].copy()
        xvals = pd.to_numeric(d[dist_col], errors="coerce").to_numpy(float)

        if dist_clip_q is not None:
            cap = float(np.nanquantile(xvals, float(dist_clip_q)))
            xvals = np.clip(xvals, a_min=None, a_max=cap)

        d["_x"] = xvals
    else:
        raise ValueError("x_kind must be 'distance' or 'steps_from_end'")

    # gaussian kernel
    def K(z):
        return np.exp(-0.5 * z * z)

    def _smooth_one_tag(resp_SxN, x_N, xg_G, h):
        S, N = resp_SxN.shape
        G = xg_G.size
        num = np.zeros((S, G), float)
        den = np.zeros((G,), float)

        for start in range(0, N, int(chunk_size)):
            stop = min(N, start + int(chunk_size))
            xs = x_N[start:stop]                       # (n,)
            rs = resp_SxN[:, start:stop]               # (S,n)
            z = (xs[:, None] - xg_G[None, :]) / h      # (n,G)
            w = K(z)                                   # (n,G)
            den += w.sum(axis=0)
            num += rs @ w

        den = den + eps
        return num / den[None, :], den

    out = {}
    tag_out = {} if return_tag_summaries else None

    for t in treat_levels:
        dt = d.loc[d[treat_col] == t].copy()
        if dt.empty:
            continue

        x_all = dt["_x"].to_numpy(float)
        xmin = float(np.nanmin(x_all))
        xmax = float(np.nanmax(x_all))
        if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
            continue

        xg = np.linspace(xmin, xmax, int(n_grid)).astype(float)
        h = float(bandwidth)
        if not np.isfinite(h) or h <= 0:
            raise ValueError("bandwidth must be positive and finite.")

        # accumulate per draw across tags, NaN-safe
        pop_sum = None                         # (S,G)
        pop_cnt = np.zeros((xg.size,), float)  # (#tags contributing at each x)
        n_tags_any = 0

        for tg, g in dt.groupby(tag_col, observed=False):
            j = g["_rowpos"].to_numpy(int)
            if j.size == 0:
                continue

            x_tg_full = g["_x"].to_numpy(float)
            x_tg = x_tg_full[np.isfinite(x_tg_full)]
            if x_tg.size < 2:
                continue

            x_min_tg = float(np.nanmin(x_tg))
            x_max_tg = float(np.nanmax(x_tg))

            resp_tg = resp[:, j]  # (S, N_tag)
            curve_SxG, den_G = _smooth_one_tag(resp_tg, x_tg_full, xg, h)

            # ---- truncation mask for this tag (same for all draws) ----
            mask = np.zeros_like(xg, dtype=bool)

            if truncate_to_data_range:
                mask |= (xg < x_min_tg) | (xg > x_max_tg)

            if truncate_by_weight:
                thr = float(den_frac) * float(np.nanmax(den_G))
                mask |= (den_G < thr)

            if mask.any():
                curve_SxG[:, mask] = np.nan

            if pop_sum is None:
                pop_sum = np.zeros_like(curve_SxG)

            good = np.isfinite(curve_SxG)  # (S,G) but mask identical across draws
            pop_sum += np.where(good, curve_SxG, 0.0)
            pop_cnt += np.isfinite(curve_SxG[0]).astype(float)

            n_tags_any += 1

            if return_tag_summaries:
                tg_mean = np.nanmean(curve_SxG, axis=0)
                tg_hdi = hdi_along_samples_2col(curve_SxG, prob=hdi_prob)
                tag_out[(t, tg)] = {
                    "x_mid": xg,
                    "mean": tg_mean,
                    "lo": tg_hdi[:, 0],
                    "hi": tg_hdi[:, 1],
                    "x_min": x_min_tg,
                    "x_max": x_max_tg,
                }

        if pop_sum is None or n_tags_any == 0:
            continue

        # equal-weight mean across tags at each x (NaN where no tags)
        denom = np.maximum(pop_cnt, 1.0)
        draw_curve = pop_sum / denom[None, :]
        draw_curve[:, pop_cnt < 1] = np.nan

        # mask unstable x where too few tags contribute (removes jagged tails)
        min_k = int(max(1, min_tags_at_x))
        unstable = pop_cnt < float(min_k)
        if unstable.any():
            draw_curve[:, unstable] = np.nan

        # ---- population post-smoothing across x ----
        if pop_post_smooth:
            bw_pop = float(pop_smooth_bandwidth) if pop_smooth_bandwidth is not None else float(bandwidth)
            w = pop_cnt.copy()
            w[unstable] = 0.0
            draw_curve = weighted_gaussian_smooth_draws(draw_curve, xg, w, bw_pop, eps=eps)
            draw_curve[:, unstable] = np.nan
            draw_curve[:, pop_cnt < 1] = np.nan

        mean = np.nanmean(draw_curve, axis=0)
        hdi = hdi_along_samples_2col(draw_curve, prob=hdi_prob)

        out[t] = {
            "label": labels.get(t, str(t)),
            "x_mid": xg,
            "mean": mean,
            "lo": hdi[:, 0],
            "hi": hdi[:, 1],
            "bandwidth": h,
            "n_tags": int(n_tags_any),
            "n_tags_by_x": pop_cnt,
            "min_tags_at_x": min_k,
            "pop_post_smooth": bool(pop_post_smooth),
            "pop_smooth_bandwidth": float(pop_smooth_bandwidth) if pop_smooth_bandwidth is not None else float(bandwidth),
        }

    if return_tag_summaries:
        return treat_levels, out, tag_out
    return treat_levels, out


def tag_curves_kernel_by_treatment(
    df,
    resp_draws,
    *,
    tag_col="tag",
    treat_col="treatment",
    treat_name_col="treatmentName",

    x_kind="distance",                 # "distance" or "steps_from_end"

    # steps-from-end x definition
    last_n_steps=150,
    steps_x_units="m",                 # "m" or "steps"
    step_length_m=60.0,

    # distance x definition
    dist_col="dist_to_dest_m",
    dist_clip_q=None,                  # None for full range; or e.g. 0.999 to cap outliers

    # smoothing controls
    n_grid=2000,
    bandwidth=240.0,
    chunk_size=4000,

    hdi_prob=0.95,
    eps=1e-12,

    # TRUNCATION CONTROLS
    truncate_to_data_range=True,
    truncate_by_weight=False,
    den_frac=0.01,
):
    """
    Returns:
      treat_levels, out_by_treat

    out_by_treat[t] = {
        "label": str,
        "x_mid": (G,),
        "bandwidth": float,
        "tags": {
            tag_value: {
                "mean": (G,), "lo": (G,), "hi": (G,),
                "n_obs": int, "x_min": float, "x_max": float
            }, ...
        }
      }
    """
    d = df.copy().reset_index(drop=True)
    n_obs = len(d)

    resp = np.asarray(resp_draws, float)
    if resp.ndim != 2 or resp.shape[1] != n_obs:
        raise ValueError(f"resp_draws must be (S, n_obs) matching df rows. Got {resp.shape} vs n_obs={n_obs}")

    d["_rowpos"] = np.arange(n_obs, dtype=int)

    if isinstance(d[treat_col].dtype, pd.CategoricalDtype):
        treat_levels = list(d[treat_col].cat.categories)
    else:
        treat_levels = sorted(d[treat_col].dropna().unique().tolist())

    if treat_name_col in d.columns:
        name_map = d.groupby(treat_col, observed=False)[treat_name_col].first().to_dict()
        labels = {t: _format_treatment_name(name_map.get(t, t)) for t in treat_levels}
    else:
        labels = {t: _format_treatment_name(t) for t in treat_levels}

    # ---- define x for each row ----
    if x_kind == "steps_from_end":
        d["_step_fwd"] = d.groupby(tag_col, observed=False).cumcount()
        d["_last_fwd"] = d.groupby(tag_col, observed=False)["_step_fwd"].transform("max")
        d["_x_steps"] = d["_step_fwd"] - d["_last_fwd"]  # last=0; earlier negative

        d = d.loc[(d["_x_steps"] >= -int(last_n_steps)) & (d["_x_steps"] <= 0)].copy()

        if str(steps_x_units).lower() in ("m", "meter", "meters"):
            d["_x"] = d["_x_steps"].astype(float) * float(step_length_m)
        elif str(steps_x_units).lower() in ("step", "steps"):
            d["_x"] = d["_x_steps"].astype(float)
        else:
            raise ValueError("steps_x_units must be 'm' or 'steps'")

    elif x_kind == "distance":
        if dist_col not in d.columns:
            raise KeyError(f"Missing distance column '{dist_col}' in df.")

        xvals = pd.to_numeric(d[dist_col], errors="coerce").to_numpy(float)
        ok = np.isfinite(xvals)
        d = d.loc[ok].copy()
        xvals = pd.to_numeric(d[dist_col], errors="coerce").to_numpy(float)

        if dist_clip_q is not None:
            cap = float(np.nanquantile(xvals, float(dist_clip_q)))
            xvals = np.clip(xvals, a_min=None, a_max=cap)

        d["_x"] = xvals
    else:
        raise ValueError("x_kind must be 'distance' or 'steps_from_end'")

    def K(z):
        return np.exp(-0.5 * z * z)

    def _smooth_one_tag(resp_SxN, x_N, xg_G, h):
        """
        resp_SxN: (S,N), x_N: (N,), xg_G: (G,)
        returns:
          curve_SxG: (S,G)
          den_G: (G,) kernel weight sum at each grid point (same for all draws)
        """
        S, N = resp_SxN.shape
        G = xg_G.size
        num = np.zeros((S, G), float)
        den = np.zeros((G,), float)

        for start in range(0, N, int(chunk_size)):
            stop = min(N, start + int(chunk_size))
            xs = x_N[start:stop]                       # (n,)
            rs = resp_SxN[:, start:stop]               # (S,n)
            z = (xs[:, None] - xg_G[None, :]) / h      # (n,G)
            w = K(z)                                   # (n,G)
            den += w.sum(axis=0)
            num += rs @ w

        den = den + eps
        return num / den[None, :], den

    out = {}

    for t in treat_levels:
        dt = d.loc[d[treat_col] == t].copy()
        if dt.empty:
            continue

        # common x-grid per treatment (align all tags in panel)
        x_all = dt["_x"].to_numpy(float)
        xmin = float(np.nanmin(x_all))
        xmax = float(np.nanmax(x_all))
        if not np.isfinite(xmin) or not np.isfinite(xmax) or xmax <= xmin:
            continue

        xg = np.linspace(xmin, xmax, int(n_grid)).astype(float)
        h = float(bandwidth)
        if not np.isfinite(h) or h <= 0:
            raise ValueError("bandwidth must be positive and finite.")

        out[t] = {
            "label": labels.get(t, str(t)),
            "x_mid": xg,
            "bandwidth": h,
            "tags": {},
        }

        for tg, g in dt.groupby(tag_col, observed=False):
            j = g["_rowpos"].to_numpy(int)
            if j.size == 0:
                continue

            x_tg = g["_x"].to_numpy(float)
            x_tg = x_tg[np.isfinite(x_tg)]
            if x_tg.size < 2:
                continue

            x_min_tg = float(np.nanmin(x_tg))
            x_max_tg = float(np.nanmax(x_tg))

            resp_tg = resp[:, j]  # (S, N_tag)
            curve_SxG, den_G = _smooth_one_tag(resp_tg, g["_x"].to_numpy(float), xg, h)

            # ---- TRUNCATION (critical): mask outside the individual's observed support ----
            mask = np.zeros_like(xg, dtype=bool)

            if truncate_to_data_range:
                mask |= (xg < x_min_tg) | (xg > x_max_tg)

            if truncate_by_weight:
                thr = float(den_frac) * float(np.nanmax(den_G))
                mask |= (den_G < thr)

            if mask.any():
                curve_SxG[:, mask] = np.nan  # breaks line & band beyond data

            tg_mean = np.nanmean(curve_SxG, axis=0)
            tg_hdi = hdi_along_samples_2col(curve_SxG, prob=hdi_prob)

            out[t]["tags"][tg] = {
                "mean": tg_mean,
                "lo": tg_hdi[:, 0],
                "hi": tg_hdi[:, 1],
                "n_obs": int(j.size),
                "x_min": x_min_tg,
                "x_max": x_max_tg,
            }

    return treat_levels, out

