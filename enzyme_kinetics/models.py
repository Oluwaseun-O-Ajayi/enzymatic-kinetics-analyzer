"""Rate laws and dose-response models.

All concentrations must be in consistent units; parameters take the units of the data.
"""

import numpy as np


def michaelis_menten(s, vmax, km):
    """Initial rate v = Vmax [S] / (Km + [S])."""
    s = np.asarray(s, dtype=float)
    return vmax * s / (km + s)


def four_parameter_logistic(x, bottom, top, log_ic50, hill):
    """Four-parameter logistic dose-response curve with the potency on a log10 scale:
    y = bottom + (top - bottom) / (1 + 10 ** (hill * (log10 x - log_ic50))).
    A positive Hill slope gives a response that falls with concentration."""
    x = np.asarray(x, dtype=float)
    with np.errstate(divide="ignore", over="ignore"):
        return bottom + (top - bottom) / (1.0 + 10.0 ** (hill * (np.log10(x) - log_ic50)))


def competitive(s, i, vmax, km, ki):
    """Competitive inhibition: v = Vmax [S] / (Km (1 + [I]/Ki) + [S])."""
    s, i = np.asarray(s, dtype=float), np.asarray(i, dtype=float)
    return vmax * s / (km * (1.0 + i / ki) + s)


def uncompetitive(s, i, vmax, km, ki):
    """Uncompetitive inhibition: v = Vmax [S] / (Km + [S] (1 + [I]/Ki))."""
    s, i = np.asarray(s, dtype=float), np.asarray(i, dtype=float)
    return vmax * s / (km + s * (1.0 + i / ki))


def mixed(s, i, vmax, km, ki, alpha):
    """Mixed inhibition: v = Vmax [S] / (Km (1 + [I]/Ki) + [S] (1 + [I]/(alpha Ki))).
    alpha = 1 gives noncompetitive inhibition."""
    s, i = np.asarray(s, dtype=float), np.asarray(i, dtype=float)
    return vmax * s / (km * (1.0 + i / ki) + s * (1.0 + i / (alpha * ki)))
