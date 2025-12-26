import numpy as np
from scipy.optimize import curve_fit

def michaelis_menten(S, Vmax, Km):
    return Vmax * S / (Km + S)

def fit_michaelis_menten(conc, velocity):
    popt, pcov = curve_fit(michaelis_menten, conc, velocity, bounds=(0, np.inf))
    Vmax, Km = popt
    return Vmax, Km, pcov  # return 3 values
