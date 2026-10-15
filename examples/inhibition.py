"""Global fit of reversible inhibition models and model comparison by AICc."""

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzyme_kinetics import compare_inhibition_models, fit_inhibition  # noqa: E402

path = ROOT / "data" / "example_inhibition.csv"
if not path.exists():
    import runpy
    runpy.run_path(str(ROOT / "examples" / "make_example_data.py"))
data = pd.read_csv(path)

ranking = compare_inhibition_models(data["substrate"], data["inhibitor"], data["velocity"])
print("Model comparison (lower AICc is better; delta below about 2 is not decisive):")
print(ranking.round(4).to_string(index=False))

best = ranking.loc[0, "model"]
fit = fit_inhibition(data["substrate"], data["inhibitor"], data["velocity"], best)
print(f"\nBest model: {best} (simulated truth: competitive, Vmax 12, Km 4, Ki 2.5)")
print(fit.table().round(4).to_string(index=False))
