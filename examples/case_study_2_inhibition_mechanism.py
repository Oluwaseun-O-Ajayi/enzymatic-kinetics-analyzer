"""Case study 2: identify an inhibition mechanism.

Question: given rates at several substrate and inhibitor concentrations, does
model comparison pick the mechanism that generated the data, and is Ki recovered?
Input:  data/example_inhibition.csv (simulated competitive inhibition,
        Vmax 12, Km 4, Ki 2.5, 4% noise, duplicates)
Output: examples/case_study_results/case2_model_comparison.csv
        examples/case_study_results/case2_competitive_fit.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS, TRUE_INHIBITION
from enzyme_kinetics import compare_inhibition_models, fit_inhibition

data = pd.read_csv(DATA / "example_inhibition.csv")
s, i, v = data["substrate"], data["inhibitor"], data["velocity"]

comparison = compare_inhibition_models(s, i, v)
comparison.to_csv(RESULTS / "case2_model_comparison.csv", index=False)
print(comparison[["model", "aicc", "delta_aicc", "akaike_weight"]].to_string(index=False))

best = comparison.iloc[0]["model"]
fit = fit_inhibition(s, i, v, model=best).table()
fit["true"] = fit["parameter"].map(TRUE_INHIBITION)
fit.to_csv(RESULTS / "case2_competitive_fit.csv", index=False)
print(f"\nSelected model: {best} (true mechanism: {TRUE_INHIBITION['mechanism']})")
print(fit.round(3).to_string(index=False))
