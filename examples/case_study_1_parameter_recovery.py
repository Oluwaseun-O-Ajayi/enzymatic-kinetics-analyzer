"""Case study 1: recover known Michaelis-Menten parameters.

Question: when the true Vmax and Km are known, how close are the fitted values,
and do the 95% confidence intervals contain the truth?
Input:  data/batch_kinetics.csv (simulated, 4% noise, eight concentrations)
Output: examples/case_study_results/case1_parameter_recovery.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, TRUE_MM
from enzyme_kinetics import batch_michaelis_menten

data = pd.read_csv(DATA / "batch_kinetics.csv")
fits = batch_michaelis_menten(data, group_by=["enzyme", "condition"],
                              substrate="substrate_concentration_uM", velocity="initial_velocity")

rows = []
for _, f in fits.iterrows():
    true_vmax, true_km = TRUE_MM[(f["enzyme"], f["condition"])]
    for p, truth in [("Vmax", true_vmax), ("Km", true_km)]:
        rows.append({"enzyme": f["enzyme"], "condition": f["condition"], "parameter": p,
                     "true": truth, "fitted": round(f[p], 3),
                     "ci_low": round(f[f"{p}_ci_low"], 3), "ci_high": round(f[f"{p}_ci_high"], 3),
                     "truth_in_ci": bool(f[f"{p}_ci_low"] <= truth <= f[f"{p}_ci_high"])})

table = pd.DataFrame(rows)
table.to_csv(RESULTS / "case1_parameter_recovery.csv", index=False)
print(table.to_string(index=False))
print(f"\nTrue value inside the 95% CI: {int(table['truth_in_ci'].sum())} of {len(table)}")
