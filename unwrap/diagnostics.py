"""ArviZ-compatible convergence diagnostics for angular posterior draws.

``rhat`` delegates existing methods to :func:`arviz.rhat`. The additional
``method="circular"`` computes a split multivariate within/between-chain
diagnostic on unit vectors. It is the square root of the Brooks--Gelman
multivariate variance ratio, not rank-normalized R-hat. See
``docs/diagnostics.md`` for its derivation, conventions, and limitations.
``circ_var_names`` selects angular variables in mixed R-hat and summary calls;
``summary`` also reports separate sine/cosine convergence and precision metrics.
"""

import logging
import warnings

import arviz as az
import numpy as np
import pandas as pd
import xarray as xr

__all__ = ["rhat", "rhat_components", "summary"]

_log = logging.getLogger("arviz")


def _same_direction(values, reference):
    """Compare constants modulo turns without erasing small within-period variation."""
    with np.errstate(over="ignore", invalid="ignore"):
        difference = values - reference
        turns = np.rint(difference / (2 * np.pi))
        residual = np.arctan2(np.sin(difference), np.cos(difference))
        tolerance = 8 * np.finfo(float).eps * np.abs(difference)
        # The tolerance applies only when removing a nonzero number of turns.
        # Its scale is the removed offset, not the absolute angular origin.
        return (difference == 0) | (
            np.isfinite(difference)
            & (turns != 0)
            & (tolerance < np.pi)
            & (np.abs(residual) <= tolerance)
        )


def _rhat_circular(ary):
    """Return the split, standard-deviation-scale multivariate PSRF for angles."""
    ary = np.asarray(ary, dtype=float)
    if not np.isfinite(ary).all():
        _log.warning("Circular R-hat requires finite angles; returning NaN.")
        return np.nan
    if ary.shape[0] < 2 or ary.shape[1] < 4:
        _log.warning(
            "Shape validation failed: input_shape: %s, minimum_shape: (chains=2, draws=4)",
            ary.shape,
        )
        return np.nan

    half = ary.shape[1] // 2
    ary = np.concatenate((ary[:, :half], ary[:, -half:]), axis=0)
    chains, draws = ary.shape

    # Circular constants may be stored with a different turn count at every draw.
    if np.all(_same_direction(ary, ary[:, :1])):
        if np.all(_same_direction(ary[:, 0], ary[0, 0])):
            return np.nan
        return np.inf

    # One common numerical reference, never a separate reference for each chain.
    with np.errstate(over="ignore"):
        phase = ary - ary[0, 0]
    if not np.isfinite(phase).all():
        reduced = np.remainder(ary, 2 * np.pi)
        phase = reduced - reduced[0, 0]
    phase = np.arctan2(np.sin(phase), np.cos(phase))

    # cos(phase)-1 has the same covariances as cos(phase), but this identity
    # preserves the small radial variation of highly concentrated directions.
    unit = np.stack((-2 * np.sin(phase / 2) ** 2, np.sin(phase)), axis=-1)
    means = unit.mean(axis=1)
    residuals = unit - means[:, None, :]
    within = np.einsum("mni,mnj->ij", residuals, residuals) / (chains * (draws - 1))
    centered_means = means - means.mean(axis=0)
    between = centered_means.T @ centered_means / (chains - 1)

    # An invertible change of coordinate scale does not change the generalized
    # eigenvalues. Balance the two coordinates before checking rank/whitening.
    scales = np.sqrt(np.diag(within))
    if np.any(scales == 0):
        _log.warning("Circular R-hat has singular within-chain covariance; returning NaN.")
        return np.nan
    scale_matrix = np.outer(scales, scales)
    within = within / scale_matrix
    between = between / scale_matrix
    if np.linalg.eigvalsh(within)[0] <= 64 * np.finfo(float).eps:
        _log.warning("Circular R-hat has singular within-chain covariance; returning NaN.")
        return np.nan

    # L^-1 S L^-T is symmetric and has the generalized eigenvalues of (S, W).
    lower = np.linalg.cholesky(within)
    whitened = np.linalg.solve(lower, between)
    whitened = np.linalg.solve(lower, whitened.T).T
    largest = np.linalg.eigvalsh((whitened + whitened.T) / 2)[-1]
    # S is positive semidefinite; remove only a possible negative roundoff error.
    largest = max(0.0, largest)
    return np.sqrt((draws - 1) / draws + (chains + 1) / chains * largest)


def _apply_diagnostic(function, data, var_names, dask_kwargs):
    """Reuse ArviZ's public conversion, selection, and labeled-array machinery."""
    if isinstance(data, np.ndarray):
        data = np.atleast_2d(data)
        if data.ndim > 2:
            raise TypeError(
                "Only uni-dimensional ndarray variables are supported. "
                "Please transform first to dataset with `az.convert_to_dataset`."
            )
        return function(data)

    dataset = az.extract(data, combined=False, var_names=var_names, keep_dataset=True)
    return az.wrap_xarray_ufunc(
        function,
        dataset,
        ufunc_kwargs={"ravel": False},
        dask_kwargs=dask_kwargs,
    )


def _circular_names(dataset, names):
    """Validate exact circular names before intersecting with a variable selection."""
    names = [names] if isinstance(names, str) else list(names)
    missing = set(names) - set(dataset.data_vars)
    if missing:
        raise KeyError(f"Circular variables not found in data: {sorted(missing)}")
    return list(dict.fromkeys(names))


def _circular_function(method):
    if method == "circular":
        return _rhat_circular
    if method in ("sin", "cos"):
        projection = np.sin if method == "sin" else np.cos
        return lambda ary: _rank_projection(ary, projection)
    raise ValueError("circ_method must be 'circular', 'sin', or 'cos'.")


def rhat(
    data, *, var_names=None, method="rank", dask_kwargs=None,
    circ_var_names=None, circ_method="circular",
):
    """Compute ArviZ R-hat or an opt-in joint circular convergence diagnostic.

    Parameters
    ----------
    data : array-like or object convertible to arviz.InferenceData
        Posterior draws. A NumPy array must have shape ``(chain, draw)``.
        Labeled inputs use the ``chain`` and ``draw`` dimensions. With
        ``method="circular"``, every selected variable must contain angles in
        radians; angular variables are not inferred automatically.
    var_names : list of str, optional
        Variables to include, using ArviZ selection rules (including ``~name``
        exclusions). Ignored for NumPy arrays, as in ``arviz.rhat``.
    method : str, default "rank"
        Existing ArviZ methods, including ``"rank"``, ``"split"``, ``"folded"``,
        ``"z_scale"``, and ``"identity"``, are passed through unchanged.
        ``"circular"`` splits each chain, embeds directions as joint unit
        vectors, and compares between- and within-chain covariance using the
        largest generalized eigenvalue. It returns the square root of the
        Brooks--Gelman multivariate variance ratio. This experimental circular
        application is not rank normalized and has no separately calibrated
        universal cutoff.
    dask_kwargs : dict, optional
        Keyword arguments for ArviZ's xarray ufunc wrapper, as in ``arviz.rhat``.
    circ_var_names : str or list of str, optional
        Exact names of angular variables in a mixed dataset. These use
        ``circ_method``; every other selected variable uses ``method``.
        Names are validated against the full posterior, then intersected with
        ``var_names``. Requires named input, not a NumPy array.
    circ_method : {"circular", "sin", "cos"}, default "circular"
        Diagnostic for ``circ_var_names``: the joint covariance statistic or
        ArviZ rank R-hat of sine or cosine. To inspect both projections, use
        ``rhat_components`` or ``summary``. Only used with ``circ_var_names``.

    Returns
    -------
    float or xarray.Dataset
        NumPy input returns a scalar; labeled/convertible inputs return a
        Dataset with sampling dimensions removed and parameter coordinates
        retained. Circular R-hat needs at least two original chains and four
        draws. Nonfinite angles, one common constant direction, or singular
        within-chain covariance yield NaN. Distinct fixed directions across
        otherwise constant split chains yield infinity.

    Notes
    -----
    The circular statistic is invariant to common rotations and whole turns,
    up to floating-point precision. It measures disagreement in mean unit
    vectors, including their lengths, and can miss different multimodal angular
    distributions with the same mean vector. It is a diagnostic, not a proof
    of convergence. Values slightly below one are retained.

    Examples
    --------
    >>> from unwrap.diagnostics import rhat
    >>> circular = rhat(idata, var_names=["mu"], method="circular")  # doctest: +SKIP
    >>> ordinary = rhat(idata, var_names=["kappa"])  # doctest: +SKIP
    >>> mixed = rhat(idata, circ_var_names=["mu", "offset"])  # doctest: +SKIP
    """
    if circ_var_names is not None and len(circ_var_names) > 0:
        if isinstance(data, np.ndarray):
            raise ValueError("circ_var_names requires named variables, not a NumPy array.")
        if method == "circular":
            raise ValueError(
                "With circ_var_names, method selects the ordinary diagnostic; "
                "use circ_method='circular' for angular variables."
            )
        function = _circular_function(circ_method)
        dataset = az.convert_to_dataset(data)
        names = _circular_names(dataset, circ_var_names)
        selected = az.extract(
            dataset, combined=False, var_names=var_names, keep_dataset=True
        )
        circular = [name for name in selected if name in names]
        ordinary = [name for name in selected if name not in names]
        parts = []
        if ordinary:
            parts.append(az.rhat(selected[ordinary], method=method, dask_kwargs=dask_kwargs))
        if circular:
            parts.append(_apply_diagnostic(function, selected[circular], None, dask_kwargs))
        if not parts:
            return az.rhat(selected, method=method, dask_kwargs=dask_kwargs)
        return xr.merge(parts)[list(selected.data_vars)]

    if method != "circular":
        return az.rhat(data, var_names=var_names, method=method, dask_kwargs=dask_kwargs)

    return _apply_diagnostic(_rhat_circular, data, var_names, dask_kwargs)


def _rank_projection(ary, projection):
    ary = np.asarray(ary)
    if not np.isfinite(ary).all():
        _log.warning("Circular R-hat requires finite angles; returning NaN.")
        return np.nan
    return az.rhat(projection(ary), method="rank")


def rhat_components(data, *, var_names=None, dask_kwargs=None):
    """Report separate ArviZ rank R-hats for sine and cosine of angular draws.

    Input conversion, variable selection, and Dask options follow ``arviz.rhat``.
    Selected variables must be angles in radians. Returns a dictionary with
    ``"sin"`` and ``"cos"`` keys; each value is an ArviZ-style scalar for NumPy
    input or Dataset for labeled/convertible input. Each diagnostic is exactly
    ArviZ's ``method="rank"`` applied to that projection. Nonfinite angles yield
    NaN, and constant projections retain ArviZ's undefined-result behavior.

    Reporting both avoids averaging away a large value in one projection.
    Each projection is periodic and smooth across the angular cut. Separate
    rank diagnostics are not jointly rotation invariant and may miss differing
    joint angular distributions with identical sine and cosine marginals.
    See ``rhat(..., method="circular")`` for the joint covariance alternative.
    """
    return {
        name: _apply_diagnostic(
            lambda ary, projection=projection: _rank_projection(ary, projection),
            data,
            var_names,
            dask_kwargs,
        )
        for name, projection in (("sin", np.sin), ("cos", np.cos))
    }


def _summary_dataset(data, group):
    """Follow ArviZ's summary group preference, including prior-only input."""
    if not isinstance(data, az.InferenceData):
        return az.convert_to_dataset(data)
    groups = data.groups()
    if group is None:
        if not groups:
            raise TypeError("InferenceData does not contain any groups")
        group = next((name for name in ("posterior", "prior") if name in groups), None)
        if group is None:
            group = groups[0]
            warnings.warn(f"Selecting first found group: {group}", stacklevel=3)
    if group not in groups:
        raise TypeError(f"InferenceData does not contain group: {group}")
    return data[group]


def _format_summary(dataset, fmt, round_to, labeller):
    """Format computed metrics using ArviZ's public label and iteration helpers."""
    fmt = fmt.lower()
    if fmt == "xarray":
        result = dataset
    elif fmt == "long":
        result = dataset.to_dataframe().reset_index().set_index("metric")
        result.index = list(result.index)
    else:
        labeller = az.labels.BaseLabeller() if labeller is None else labeller
        rows, labels = [], []
        for name, selection, indices, values in az.sel_utils.xarray_var_iter(
            dataset, skip_dims={"metric"}
        ):
            rows.append(values)
            labels.append(labeller.make_label_flat(name, selection, indices))
        result = pd.DataFrame(rows, index=labels, columns=dataset.metric.values)
    if round_to is not None and round_to not in ("none", "None"):
        return result.round(round_to)
    if round_to is None and fmt != "xarray":
        decimals = {}
        for column in result.columns:
            # Projection diagnostics have the same precision as their ArviZ counterparts.
            metric = column
            if fmt == "wide" and column.endswith(("_sin", "_cos")):
                metric = column[:-4]
            decimals[column] = 2 if metric == "r_hat" else 0 if metric in (
                "ess_bulk", "ess_tail"
            ) else 3
        return result.round(decimals)
    return result


def summary(
    data, var_names=None, filter_vars=None, group=None, fmt="wide", kind="all",
    round_to=None, circ_var_names=None, stat_focus="mean", stat_funcs=None,
    extend=True, hdi_prob=None, skipna=False, labeller=None, coords=None,
    index_origin=None, order=None, *, circ_method="circular",
    include_circular_diagnostics=True,
):
    """Summarize mixed circular and ordinary variables with ArviZ-style output.

    Accepts the arguments of :func:`arviz.summary`. With no circular variables,
    delegates to ArviZ unchanged. ``circ_var_names`` contains exact variable
    names, in radians; names excluded by ``var_names`` are ignored. Circular
    summaries require ``stat_focus="mean"``, as in ArviZ.

    For circular rows, ``r_hat`` uses ``circ_method`` (``"circular"``, ``"sin"``,
    or ``"cos"``; see :func:`rhat`). Mean, SD, HDI and circular mean MCSE come
    from ArviZ. The raw-angle ``ess_bulk``, ``ess_tail`` and ``mcse_sd`` cells
    are NaN because these are not circular diagnostics. In particular, no
    universal circular tail ESS or MCSE of circular SD is claimed here.

    When ``include_circular_diagnostics=True`` (default), append all five ArviZ
    diagnostic metrics for each sine/cosine projection, with ``_sin``/``_cos``
    suffixes. These cells are NaN for ordinary rows. Projected MCSEs describe
    the mean/SD of the projection, not an angular mean/SD. Both projection
    R-hats are reported separately, never averaged. ``kind="stats"`` and
    ``extend=False`` do not add diagnostics. Ordinary rows retain ArviZ's
    statistics, and all three output formats and custom labels are supported.

    The default joint circular R-hat is experimental, not rank-normalized;
    its threshold is not calibrated to the usual rank R-hat cutoff. Neither
    it nor the projected diagnostics can prove convergence. Also inspect
    trace plots, NUTS divergences, energy and effective sample sizes.

    Examples
    --------
    >>> from unwrap.diagnostics import summary
    >>> summary(idata, circ_var_names=["Intercept", "angular_offset"])  # doctest: +SKIP
    >>> summary(idata, circ_var_names=["mu"], kind="diagnostics")  # doctest: +SKIP
    """
    kwargs = dict(
        var_names=var_names, filter_vars=filter_vars, group=group, fmt=fmt,
        kind=kind, round_to=round_to, circ_var_names=circ_var_names,
        stat_focus=stat_focus, stat_funcs=stat_funcs, extend=extend,
        hdi_prob=hdi_prob, skipna=skipna, labeller=labeller, coords=coords,
        index_origin=index_origin, order=order,
    )
    if circ_var_names is None or len(circ_var_names) == 0:
        return az.summary(data, **kwargs)

    dataset = _summary_dataset(data, group)
    names = _circular_names(dataset, circ_var_names)
    selected = az.extract(
        dataset, combined=False, var_names=var_names, filter_vars=filter_vars,
        keep_dataset=True,
    )
    selected = az.utils.get_coords(selected, {} if coords is None else coords)
    names = [name for name in selected if name in names]
    # Leave validation and all non-circular summary behavior with ArviZ.
    if not names:
        kwargs["circ_var_names"] = None
        return az.summary(data, **kwargs)
    if kind == "stats" or not extend:
        kwargs.update(var_names=None, filter_vars=None, group=None, coords=None,
                      circ_var_names=names)
        return az.summary(selected, **kwargs)

    function = _circular_function(circ_method)
    # Validate the requested format even though metrics are computed as xarray.
    if not isinstance(fmt, str) or fmt.lower() not in ("wide", "long", "xarray"):
        raise TypeError(f"Invalid format: '{fmt}'. Formatting options are: ('wide', 'long', 'xarray')")
    kwargs.update(var_names=None, filter_vars=None, group=None, coords=None,
                  circ_var_names=names, fmt="xarray", round_to="none")
    result = az.summary(selected, **kwargs)
    if len(set(result.metric.values)) != len(result.metric):
        raise ValueError("Custom statistic names must not duplicate summary metric names.")
    circ_data = selected[names]
    circular_rhat = _apply_diagnostic(function, circ_data, None, None)
    # ArviZ 0.x leaves mcse_mean linear for kind='diagnostics'; consistently
    # obtain its circular mean MC error from the mean-focused all summary.
    circular_mcse = result.sel(metric="mcse_mean")
    if kind == "diagnostics":
        circular_mcse = az.summary(
            circ_data, circ_var_names=names, fmt="xarray", round_to="none",
            hdi_prob=hdi_prob, skipna=skipna,
        ).sel(metric="mcse_mean")
    for name in names:
        result[name].loc[dict(metric="r_hat")] = circular_rhat[name]
        result[name].loc[dict(metric="mcse_mean")] = circular_mcse[name]
        result[name].loc[dict(metric=["ess_bulk", "ess_tail", "mcse_sd"])] = np.nan

    if include_circular_diagnostics:
        extras = []
        for suffix, projection in (("sin", np.sin), ("cos", np.cos)):
            projected = az.summary(
                projection(circ_data), kind="diagnostics", fmt="xarray", round_to="none"
            )
            metrics = [f"{metric}_{suffix}" for metric in projected.metric.values]
            if set(metrics) & set(result.metric.values):
                raise ValueError("Custom statistic names collide with circular diagnostic names.")
            projected = projected.assign_coords(metric=metrics)
            padded = xr.Dataset({
                name: xr.full_like(result[name].isel(metric=0, drop=True), np.nan)
                .expand_dims(metric=metrics)
                for name in result.data_vars
            })
            for name in names:
                padded[name] = projected[name]
            extras.append(padded)
        result = xr.concat([result, *extras], dim="metric")
    return _format_summary(result, fmt, round_to, labeller)
