import numpy as np


def parameter_confidence_intervals(popt, pcov, confidence=0.95):
    z = 1.96  # 95% confidence interval
    intervals = []

    for i, param in enumerate(popt):
        sigma = np.sqrt(pcov[i, i])
        intervals.append((param - z * sigma, param + z * sigma))

    return intervals
