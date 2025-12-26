import numpy as np
from scipy.optimize import curve_fit


def michaelis_menten(substrate_conc, vmax, km):
    return (vmax * substrate_conc) / (km + substrate_conc)


def fit_michaelis_menten(substrate_conc, velocity):
    popt, pcov = curve_fit(
        michaelis_menten,
        substrate_conc,
        velocity,
        p0=[np.max(velocity), np.median(substrate_conc)]
    )
    return popt, pcov
