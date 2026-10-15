"""Fit the Michaelis-Menten equation to one data set and to a batch of enzymes and conditions."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzyme_kinetics import batch_michaelis_menten, fit_michaelis_menten  # noqa: E402
from enzyme_kinetics.plotting import plot_michaelis_menten  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

data = pd.read_csv(ROOT / "data" / "example_michaelis_menten.csv")
fit = fit_michaelis_menten(data["substrate"], data["velocity"])
print("Single data set (95% confidence intervals, t distribution):")
print(fit.table().round(4).to_string(index=False))
ax = plot_michaelis_menten(fit, data["substrate"], data["velocity"])
ax.figure.tight_layout()
ax.figure.savefig(OUT / "michaelis_menten.png", dpi=200)
plt.close(ax.figure)

batch = pd.read_csv(ROOT / "data" / "batch_kinetics.csv")
table = batch_michaelis_menten(batch, ["enzyme", "condition"],
                               substrate="substrate_concentration_uM", velocity="initial_velocity")
print("\nBatch fit:")
print(table.round(4).to_string(index=False))
table.to_csv(OUT / "batch_kinetics_results.csv", index=False)
