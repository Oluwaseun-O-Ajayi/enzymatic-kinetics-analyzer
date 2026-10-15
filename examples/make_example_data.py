"""Create the simulated example data sets (fixed seed, known parameters):

- data/batch_kinetics.csv: two enzymes x two conditions, eight substrate concentrations,
  4% multiplicative noise. True (Vmax, Km): EnzA/pH7 (1.0, 5.0), EnzA/pH8 (1.4, 3.0),
  EnzB/WT (0.8, 10.0), EnzB/Mutant (0.8, 25.0).
- data/example_inhibition.csv: competitive inhibition with Vmax 12, Km 4, Ki 2.5,
  five inhibitor and eight substrate concentrations in duplicate, 4% noise.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzyme_kinetics import competitive, michaelis_menten  # noqa: E402

rng = np.random.default_rng(13)
S = (0.5, 1, 2, 4, 8, 16, 32, 64)

rows = []
truth = {("EnzA", "pH7"): (1.0, 5.0), ("EnzA", "pH8"): (1.4, 3.0),
         ("EnzB", "WT"): (0.8, 10.0), ("EnzB", "Mutant"): (0.8, 25.0)}
for (enzyme, condition), (vmax, km) in truth.items():
    for s in S:
        v = michaelis_menten(s, vmax, km) * (1 + rng.normal(0, 0.04))
        rows.append({"enzyme": enzyme, "condition": condition,
                     "substrate_concentration_uM": s, "initial_velocity": round(float(v), 4)})
pd.DataFrame(rows).to_csv(ROOT / "data" / "batch_kinetics.csv", index=False)

rows = []
for i in (0.0, 1.0, 2.5, 5.0, 10.0):
    for s in S:
        for rep in (1, 2):
            v = competitive(s, i, 12.0, 4.0, 2.5) * (1 + rng.normal(0, 0.04))
            rows.append({"substrate": s, "inhibitor": i, "replicate": rep, "velocity": round(float(v), 4)})
pd.DataFrame(rows).to_csv(ROOT / "data" / "example_inhibition.csv", index=False)
print("Wrote data/batch_kinetics.csv and data/example_inhibition.csv")
