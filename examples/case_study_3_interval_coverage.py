"""Case study 3: why the confidence intervals use the t distribution.

Question: with six or eight substrate concentrations, how often does a nominal
95% interval for Km contain the true value, using Student's t (as the package
does) versus the normal approximation (z = 1.96)?
Input:  simulated experiments (Vmax 1, Km 5, 4% multiplicative noise, fixed seed)
Output: examples/case_study_results/case3_interval_coverage.csv
"""

import numpy as np
import pandas as pd
from _case_study_common import RESULTS
from enzyme_kinetics import fit_michaelis_menten, michaelis_menten

VMAX, KM, N_SIM = 1.0, 5.0, 1000
designs = {6: np.array([0.5, 1, 4, 8, 32, 64]), 8: np.array([0.5, 1, 2, 4, 8, 16, 32, 64])}
rng = np.random.default_rng(2026)

rows = []
for n_points, s in designs.items():
    t_hits = z_hits = 0
    for _ in range(N_SIM):
        v = michaelis_menten(s, VMAX, KM) * (1 + rng.normal(0, 0.04, len(s)))
        fit = fit_michaelis_menten(s, v)
        t_hits += fit.ci_low["Km"] <= KM <= fit.ci_high["Km"]
        z_hits += abs(fit.params["Km"] - KM) <= 1.96 * fit.stderr["Km"]
    rows.append({"points": n_points, "simulations": N_SIM,
                 "coverage_t": t_hits / N_SIM, "coverage_normal": z_hits / N_SIM, "nominal": 0.95})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case3_interval_coverage.csv", index=False)
print(table.to_string(index=False))
