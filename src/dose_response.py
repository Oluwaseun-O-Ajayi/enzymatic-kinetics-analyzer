import numpy as np
from scipy.optimize import curve_fit

# 4-parameter logistic (4PL) function
def logistic4(x, bottom, top, ic50, hill_slope):
    return bottom + (top - bottom) / (1 + (x/ic50)**hill_slope)

# IC50 fitting function
def fit_ic50(x, y, return_curve=False):
    """
    Fit 4PL dose-response curve.
    
    Parameters
    ----------
    x : array-like
        Inhibitor concentrations
    y : array-like
        Response values
    return_curve : bool
        If True, returns fitted curve instead of parameters
    
    Returns
    -------
    popt : array
        Fitted parameters [bottom, top, IC50, hill_slope] if return_curve=False
    y_fit : array
        Fitted curve if return_curve=True
    """
    # Fit 4PL curve
    popt, _ = curve_fit(logistic4, x, y, bounds=(0, [100, 100, np.inf, 5]))
    
    if return_curve:
        return logistic4(x, *popt)
    else:
        return popt
