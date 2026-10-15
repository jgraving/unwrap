# Circular and positive links

In a circular GLM or GLMM, a linear predictor combines an intercept with
effects of treatment, density, or other predictors. An **inverse link** turns
that sum into a parameter of the von Mises distribution. For example,
`circular_modulo` turns a predicted angle into a direction on the circle;
`softplus` turns a concentration predictor into a positive `kappa`. This lets
the same regression framework describe both reorientation and changes in
directional consistency. See the [model guide](models.md) for GLM and GLMM
formulas and their biological interpretation.

In the API, `link` maps a response parameter to predictor scale; `linkinv`
maps the linear predictor to response-parameter scale. Softplus is an inverse
link; inverse softplus reverses that mapping for positive inputs. A
concentration coefficient is therefore interpreted before softplus, rather
than as a constant additive change in `kappa`.

## Backends

| Import | Provides |
| --- | --- |
| `unwrap.link` or `unwrap.link.numpy` | NumPy functions for scalars and arrays |
| `unwrap.link.pymc` | Functions usable in symbolic PyTensor expressions |
| `unwrap.link.bambi.LINKS` | Bambi Link objects combining both backends |

`NUMPY_LINKS`, `PYMC_LINKS`, and `BAMBI_LINKS` are also available from
`unwrap.link`. Importing that package currently loads all three backends.

```python
import numpy as np
from unwrap.link import circular_modulo
from unwrap.link.bambi import LINKS

angles = circular_modulo(np.deg2rad(np.array([179, 181])))
links = {"mu": LINKS["circular_modulo"], "kappa": LINKS["softplus"]}
# Pass links to bmb.Model(..., family="vonmises", link=links).
```

## Names and aliases

| Registry name | Inverse-link operation | Aliases |
| --- | --- | --- |
| `softplus` | `log(1 + exp(eta))` | `logexpm1`, `log1pexp` |
| `reciprocal_softplus` | `1 / softplus(eta)` | `inv_softplus` |
| `modulo` | Wrap to `[0, 2*pi)` | `mod` |
| `circular_modulo` | Wrap to `[-pi, pi)` | `circ_mod`, `shifted_modulo`, `shift_mod`, `rotated_modulo`, `rot_mod` |
| `rectifier` | Floor at float64 machine epsilon | `relu` |
| `atan2` | `atan2(sin(eta), cos(eta))` | None |

These are registry keys. Not every alias is exported as a standalone NumPy
function. Bambi exports Link objects under every registry name.

`inv_softplus` names the **reciprocal** softplus link. The function
`inverse_softplus` is the mathematical inverse of softplus. Their meanings differ.
`shifted_modulo` and `rotated_modulo` are exact aliases of `circular_modulo`,
without any adjustable shift or rotation argument.

## Boundaries and numerical limits

Both -pi and pi map to -pi under `circular_modulo`. `atan2` may give either boundary
sign because of floating-point rounding. All circular wrappers are many-to-one,
with discontinuities at the branch cut. They have no global mathematical inverse;
the registries use the wrapper in both slots as a modelling convention.

The rectifier also has no inverse and floors at epsilon rather than zero.
Softplus and its inverse use stable formulas, but finite-precision limits still
apply: very negative inputs can underflow to zero. The inverse requires a positive
input for a finite real result, and the reciprocal can overflow near zero.
