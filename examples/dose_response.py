"""IC50 analysis with replicates, for one compound and for a batch of compounds."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzyme_kinetics import (batch_dose_response, fit_dose_response, read_compound_table,  # noqa: E402
                             read_replicate_table)
from enzyme_kinetics.plotting import plot_dose_response  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

# One compound, three replicates per concentration: every replicate enters the fit.
reps = read_replicate_table(ROOT / "data" / "example_ic50.csv")
fit = fit_dose_response(reps["concentration"], reps["response"])
print("Single compound (all replicates fitted):")
print(fit.table().round(4).to_string(index=False))
ax = plot_dose_response(fit, reps["concentration"], reps["response"])
ax.figure.tight_layout()
ax.figure.savefig(OUT / "dose_response.png", dpi=200)
plt.close(ax.figure)

for name in ("batch_ic50.csv", "batch_ic50_with_reps.csv"):
    long = read_compound_table(ROOT / "data" / name)
    table = batch_dose_response(long)
    print(f"\nBatch fit of {name}:")
    print(table[["compound", "IC50", "IC50_ci_low", "IC50_ci_high", "hill", "n", "note"]].round(4).to_string(index=False))
    table.to_csv(OUT / name.replace(".csv", "_results.csv"), index=False)

# The replicate batch stops at 0.5, before the lower plateau is reached. For percent-activity
# data, fixing the plateaus at 0 and 100 makes the IC50 estimate depend less on extrapolation.
long = read_compound_table(ROOT / "data" / "batch_ic50_with_reps.csv")
fixed = batch_dose_response(long, fix_bottom=0.0, fix_top=100.0)
print("\nSame batch with bottom fixed at 0 and top at 100:")
print(fixed[["compound", "IC50", "IC50_ci_low", "IC50_ci_high", "hill", "n", "note"]].round(4).to_string(index=False))
