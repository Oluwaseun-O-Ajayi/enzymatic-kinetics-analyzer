"""Shared helpers for the case-study scripts."""

from pathlib import Path
import sys
import warnings

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore", category=RuntimeWarning)

DATA = ROOT / "data"
RESULTS = Path(__file__).resolve().parent / "case_study_results"
RESULTS.mkdir(exist_ok=True)

# True parameters used to simulate the bundled data (see examples/make_example_data.py)
TRUE_MM = {("EnzA", "pH7"): (1.0, 5.0), ("EnzA", "pH8"): (1.4, 3.0),
           ("EnzB", "WT"): (0.8, 10.0), ("EnzB", "Mutant"): (0.8, 25.0)}
TRUE_INHIBITION = {"Vmax": 12.0, "Km": 4.0, "Ki": 2.5, "mechanism": "competitive"}
