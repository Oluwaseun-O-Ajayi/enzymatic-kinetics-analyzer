"""Nonlinear least-squares fitting with t-based confidence intervals and model comparison.

Standard errors come from the covariance matrix returned by scipy.optimize.curve_fit,
scaled by the residual variance. Confidence intervals use Student's t distribution with
n - p degrees of freedom, which matters for the small data sets typical of kinetics.
These are asymptotic (Wald) intervals; they can be too narrow when the data constrain a
parameter poorly, for example Km when no substrate concentration exceeds Km.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit

from . import models


@dataclass
class FitResult:
    """Parameters, uncertainties and goodness-of-fit statistics of one fit."""

    model: str
    params: Dict[str, float]
    stderr: Dict[str, float]
    ci_low: Dict[str, float]
    ci_high: Dict[str, float]
    n: int
    rss: float
    r_squared: float
    aic: float
    aicc: float
    confidence: float
    residuals: np.ndarray = field(repr=False)
    n_params: int = 0
    derived: Dict[str, Dict[str, float]] = field(default_factory=dict)

    @property
    def dof(self) -> int:
        """Residual degrees of freedom: points minus fitted (not fixed) parameters."""
        return self.n - self.n_params

    def table(self) -> pd.DataFrame:
        """One row per parameter (and per derived quantity) with value, SE and CI."""
        rows = [{"parameter": k, "value": v, "stderr": self.stderr[k],
                 "ci_low": self.ci_low[k], "ci_high": self.ci_high[k]} for k, v in self.params.items()]
        for k, d in self.derived.items():
            rows.append({"parameter": k, **d})
        return pd.DataFrame(rows)


def _fit(func, xdata, y, p0, bounds, names, model, confidence, sigma=None):
    y = np.asarray(y, dtype=float)
    n, p = len(y), len(p0)
    if n <= p:
        raise ValueError(f"{model}: need more than {p} data points, got {n}")
    popt, pcov = curve_fit(func, xdata, y, p0=p0, bounds=bounds, sigma=sigma,
                           absolute_sigma=False, maxfev=20000)
    resid = y - func(xdata, *popt)
    w = 1.0 if sigma is None else 1.0 / np.asarray(sigma, dtype=float) ** 2
    rss = float(np.sum(w * resid ** 2))
    tss = float(np.sum(w * (y - np.average(y, weights=np.broadcast_to(w, y.shape))) ** 2))
    se = np.sqrt(np.clip(np.diag(pcov), 0, None))
    tcrit = stats.t.ppf(0.5 + confidence / 2.0, n - p)
    k = p + 1  # parameters plus the error variance
    aic = n * np.log(max(rss, 1e-300) / n) + 2 * k
    aicc = aic + (2 * k * (k + 1)) / (n - k - 1) if n - k - 1 > 0 else np.inf
    return FitResult(
        model=model,
        params=dict(zip(names, map(float, popt))),
        stderr=dict(zip(names, map(float, se))),
        ci_low=dict(zip(names, map(float, popt - tcrit * se))),
        ci_high=dict(zip(names, map(float, popt + tcrit * se))),
        n=n, rss=rss, r_squared=1.0 - rss / tss if tss > 0 else np.nan,
        aic=float(aic), aicc=float(aicc), confidence=confidence, residuals=resid, n_params=p,
    )


def fit_michaelis_menten(substrate, velocity, sigma=None, confidence=0.95) -> FitResult:
    """Fit v = Vmax [S] / (Km + [S]) to initial rates.

    Args:
        substrate, velocity: arrays of equal length; replicate measurements may be
            passed as repeated substrate values (all points are fitted).
        sigma: optional standard deviations of the velocities (weighted fit).
        confidence: confidence level of the intervals.
    """
    s = np.asarray(substrate, dtype=float)
    v = np.asarray(velocity, dtype=float)
    p0 = [float(np.max(v)), float(np.median(s))]
    return _fit(models.michaelis_menten, s, v, p0, ([0, 0], [np.inf, np.inf]),
                ["Vmax", "Km"], "michaelis_menten", confidence, sigma)


def fit_dose_response(concentration, response, sigma=None, confidence=0.95,
                      fix_bottom: Optional[float] = None, fix_top: Optional[float] = None) -> FitResult:
    """Fit a four-parameter logistic curve and report IC50 with an asymmetric interval.

    The potency is fitted as log10(IC50), so its interval is symmetric on the log scale
    and the reported IC50 interval (10 ** limits) is asymmetric, as is usual.
    Concentrations must be positive. Bottom and/or top may be fixed (for example to 0
    and 100 for percent activity), which turns the model into a three- or two-parameter fit.
    """
    x = np.asarray(concentration, dtype=float)
    y = np.asarray(response, dtype=float)
    if np.any(x <= 0):
        raise ValueError("Dose-response concentrations must be positive")
    lo, hi = float(np.min(y)), float(np.max(y))
    span = hi - lo if hi > lo else 1.0
    guess = {"bottom": lo, "top": hi, "log_ic50": float(np.median(np.log10(x))), "hill": 1.0}
    lower = {"bottom": lo - span, "top": lo - span, "log_ic50": np.log10(x.min()) - 3, "hill": 0.05}
    upper = {"bottom": hi + span, "top": hi + span, "log_ic50": np.log10(x.max()) + 3, "hill": 10.0}
    fixed = {"bottom": fix_bottom, "top": fix_top}
    free = [k for k in ("bottom", "top", "log_ic50", "hill") if fixed.get(k) is None]

    def f(xx, *theta):
        vals = dict(zip(free, theta))
        vals.update({k: v for k, v in fixed.items() if v is not None})
        return models.four_parameter_logistic(xx, vals["bottom"], vals["top"], vals["log_ic50"], vals["hill"])

    res = _fit(f, x, y, [guess[k] for k in free], ([lower[k] for k in free], [upper[k] for k in free]),
               free, "four_parameter_logistic", confidence, sigma)
    for k, v in fixed.items():
        if v is not None:
            res.params[k] = float(v)
            res.stderr[k] = np.nan
            res.ci_low[k] = res.ci_high[k] = float(v)
    li = res.params["log_ic50"]
    res.derived["IC50"] = {"value": 10 ** li, "stderr": np.nan,
                           "ci_low": 10 ** res.ci_low["log_ic50"], "ci_high": 10 ** res.ci_high["log_ic50"]}
    return res


_INHIBITION = {
    "competitive": (models.competitive, ["Vmax", "Km", "Ki"]),
    "uncompetitive": (models.uncompetitive, ["Vmax", "Km", "Ki"]),
    "mixed": (models.mixed, ["Vmax", "Km", "Ki", "alpha"]),
}


def fit_inhibition(substrate, inhibitor, velocity, model="competitive", sigma=None, confidence=0.95) -> FitResult:
    """Global fit of a reversible inhibition model to rates measured at several substrate
    and inhibitor concentrations ("competitive", "uncompetitive" or "mixed")."""
    if model not in _INHIBITION:
        raise ValueError(f"model must be one of {sorted(_INHIBITION)}")
    s = np.asarray(substrate, dtype=float)
    i = np.asarray(inhibitor, dtype=float)
    v = np.asarray(velocity, dtype=float)
    func, names = _INHIBITION[model]
    i_pos = i[i > 0]
    p0 = {"Vmax": float(v.max()), "Km": float(np.median(s)),
          "Ki": float(np.median(i_pos)) if len(i_pos) else 1.0, "alpha": 1.0}
    lower = [0, 0, 0, 1e-3][:len(names)]
    upper = [np.inf, np.inf, np.inf, 1e3][:len(names)]
    return _fit(lambda X, *p: func(X[0], X[1], *p), np.vstack([s, i]), v,
                [p0[n] for n in names], (lower, upper), names, model, confidence, sigma)


def compare_inhibition_models(substrate, inhibitor, velocity,
                              models_to_fit: Sequence[str] = ("competitive", "uncompetitive", "mixed"),
                              sigma=None) -> pd.DataFrame:
    """Fit several inhibition models and rank them by small-sample AIC (AICc).

    Returns one row per model with its parameters, RSS, AICc, the difference from the best
    model and the Akaike weight. A difference below about 2 means the data do not
    distinguish the models."""
    rows: List[dict] = []
    for m in models_to_fit:
        r = fit_inhibition(substrate, inhibitor, velocity, m, sigma)
        rows.append({"model": m, **r.params, "rss": r.rss, "aicc": r.aicc})
    df = pd.DataFrame(rows)
    df["delta_aicc"] = df["aicc"] - df["aicc"].min()
    w = np.exp(-0.5 * df["delta_aicc"])
    df["akaike_weight"] = w / w.sum()
    return df.sort_values("aicc").reset_index(drop=True)
