"""Fit many data sets at once and collect the results in one table."""

from typing import Sequence, Union

import numpy as np
import pandas as pd

from .fitting import FitResult, fit_dose_response, fit_michaelis_menten


def _row(keys, result: FitResult):
    row = dict(keys)
    for name, d in result.table().set_index("parameter").iterrows():
        row[name] = d["value"]
        row[f"{name}_ci_low"] = d["ci_low"]
        row[f"{name}_ci_high"] = d["ci_high"]
    row["n"] = result.n
    row["r_squared"] = result.r_squared
    return row


def batch_michaelis_menten(df: pd.DataFrame, group_by: Union[str, Sequence[str]],
                           substrate: str = "substrate", velocity: str = "velocity",
                           confidence: float = 0.95) -> pd.DataFrame:
    """Fit the Michaelis-Menten equation to each group of a long-format table."""
    group_by = [group_by] if isinstance(group_by, str) else list(group_by)
    rows = []
    for keys, g in df.groupby(group_by, sort=False):
        keys = keys if isinstance(keys, tuple) else (keys,)
        if len(g) <= 2:
            rows.append({**dict(zip(group_by, keys)), "n": len(g),
                         "note": "too few points to fit two parameters"})
            continue
        res = fit_michaelis_menten(g[substrate], g[velocity], confidence=confidence)
        rows.append({**_row(zip(group_by, keys), res), "note": ""})
    return pd.DataFrame(rows)


def batch_dose_response(df: pd.DataFrame, group_by: str = "compound", concentration: str = "concentration",
                        response: str = "response", confidence: float = 0.95, **kwargs) -> pd.DataFrame:
    """Fit a four-parameter logistic curve to each group of a long-format table.
    A note is added when the fitted IC50 lies outside the tested concentration range."""
    rows = []
    for key, g in df.groupby(group_by, sort=False):
        res = fit_dose_response(g[concentration], g[response], confidence=confidence, **kwargs)
        row = _row([(group_by, key)], res)
        row["note"] = "; ".join(dose_response_warnings(res, g[concentration], g[response]))
        rows.append(row)
    return pd.DataFrame(rows)


def dose_response_warnings(result: FitResult, concentration, response) -> list:
    """Plain-language warnings for dose-response fits that the data do not support well."""
    x, y = np.asarray(concentration, float), np.asarray(response, float)
    span = y.max() - y.min() if y.max() > y.min() else 1.0
    p, out = result.params, []
    ic50 = result.derived["IC50"]["value"]
    levels = np.unique(x)
    if not x.min() <= ic50 <= x.max():
        out.append("IC50 outside tested range")
    elif (levels < ic50).sum() < 2 or (levels > ic50).sum() < 2:
        out.append("fewer than two concentrations on one side of the IC50")
    free = lambda k: not np.isnan(result.stderr.get(k, np.nan))  # fixed plateaus are not checked
    if (free("bottom") and p["bottom"] < y.min() - 0.1 * span) or (free("top") and p["top"] > y.max() + 0.1 * span):
        out.append("plateau extrapolated beyond the data; consider fixing bottom or top")
    if p["hill"] > 3:
        out.append("steep Hill slope; check the concentration spacing")
    return out


def summarize_replicates(df: pd.DataFrame, by: Sequence[str], value: str = "response") -> pd.DataFrame:
    """Mean, standard deviation and count of replicate measurements (for plotting)."""
    out = df.groupby(list(by))[value].agg(["mean", "std", "count"]).reset_index()
    out["std"] = out["std"].fillna(0.0)
    return out
