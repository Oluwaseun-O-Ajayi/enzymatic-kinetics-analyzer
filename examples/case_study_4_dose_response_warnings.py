"""Case study 4: dose-response fits that need a second look.

Question: for three compounds measured in triplicate, which IC50 estimates can
be trusted as reported, and what does the package flag for the others?
Input:  data/batch_ic50_with_reps.csv (illustrative values; concentrations stop
        at 0.5, so some curves are not fully defined)
Output: examples/case_study_results/case4_dose_response.csv
"""

import pandas as pd
from _case_study_common import DATA, RESULTS
from enzyme_kinetics import batch_dose_response, read_compound_table

long = read_compound_table(DATA / "batch_ic50_with_reps.csv")
fits = batch_dose_response(long)
cols = ["compound", "IC50", "IC50_ci_low", "IC50_ci_high", "hill", "bottom", "top", "n", "note"]
table = fits[cols].round(3)
table.to_csv(RESULTS / "case4_dose_response.csv", index=False)
print(table.to_string(index=False))
