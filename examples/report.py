"""Run the example analyses and collect tables and figures in one PDF report."""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from enzyme_kinetics import fit_dose_response, fit_michaelis_menten, read_replicate_table  # noqa: E402
from enzyme_kinetics.plotting import plot_dose_response, plot_michaelis_menten  # noqa: E402
from enzyme_kinetics.report import write_pdf_report  # noqa: E402

OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)

mm = pd.read_csv(ROOT / "data" / "example_michaelis_menten.csv")
mm_fit = fit_michaelis_menten(mm["substrate"], mm["velocity"])
ic = read_replicate_table(ROOT / "data" / "example_ic50.csv")
ic_fit = fit_dose_response(ic["concentration"], ic["response"])

path = write_pdf_report(OUT / "enzyme_kinetics_report.pdf", "Enzyme kinetics analysis", [
    ("Michaelis-Menten fit", mm_fit.table(), plot_michaelis_menten(mm_fit, mm["substrate"], mm["velocity"]).figure),
    ("Dose-response fit", ic_fit.table(), plot_dose_response(ic_fit, ic["concentration"], ic["response"]).figure),
])
print(f"Report written to {path}")
