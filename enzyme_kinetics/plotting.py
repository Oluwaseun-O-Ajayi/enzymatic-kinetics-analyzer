"""Plots of data and fitted curves."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import models
from .fitting import FitResult


def _means(x, y):
    d = pd.DataFrame({"x": np.asarray(x, float), "y": np.asarray(y, float)})
    g = d.groupby("x")["y"].agg(["mean", "std", "count"]).reset_index()
    return g["x"].to_numpy(), g["mean"].to_numpy(), g["std"].fillna(0).to_numpy(), g["count"].max()


def plot_michaelis_menten(result: FitResult, substrate, velocity, ax=None, label=None):
    """Data (mean +/- SD when replicated) and the fitted Michaelis-Menten curve."""
    ax = ax or plt.subplots(figsize=(5, 4))[1]
    xs, ym, ysd, nrep = _means(substrate, velocity)
    ax.errorbar(xs, ym, yerr=ysd if nrep > 1 else None, fmt="o", capsize=3, label=label or "Data")
    grid = np.linspace(0, xs.max() * 1.1, 200)
    p = result.params
    ax.plot(grid, models.michaelis_menten(grid, p["Vmax"], p["Km"]), "-",
            label=f"Fit: Vmax {p['Vmax']:.3g}, Km {p['Km']:.3g}")
    ax.set_xlabel("Substrate concentration")
    ax.set_ylabel("Initial rate")
    ax.legend(frameon=False, fontsize=8)
    return ax


def plot_dose_response(result: FitResult, concentration, response, ax=None, label=None):
    """Data (mean +/- SD when replicated) and the fitted four-parameter logistic curve."""
    ax = ax or plt.subplots(figsize=(5, 4))[1]
    xs, ym, ysd, nrep = _means(concentration, response)
    ax.errorbar(xs, ym, yerr=ysd if nrep > 1 else None, fmt="o", capsize=3, label=label or "Data")
    grid = np.logspace(np.log10(xs.min()) - 0.5, np.log10(xs.max()) + 0.5, 200)
    p = result.params
    ic50 = result.derived["IC50"]
    ax.plot(grid, models.four_parameter_logistic(grid, p["bottom"], p["top"], p["log_ic50"], p["hill"]), "-",
            label=f"IC50 {ic50['value']:.3g} ({ic50['ci_low']:.3g} to {ic50['ci_high']:.3g})")
    ax.set_xscale("log")
    ax.set_xlabel("Inhibitor concentration")
    ax.set_ylabel("Response")
    ax.legend(frameon=False, fontsize=8)
    return ax
