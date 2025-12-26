import numpy as np
from scipy.optimize import curve_fit


def ic50_model(inhibitor_conc, top, bottom, ic50):
    return bottom + (top - bottom) / (1 + (inhibitor_conc / ic50))


def fit_ic50(inhibitor_conc, response):
    popt, pcov = curve_fit(
        ic50_model,
        inhibitor_conc,
        response,
        p0=[np.max(response), np.min(response), np.median(inhibitor_conc)]
    )
    return popt, pcov
