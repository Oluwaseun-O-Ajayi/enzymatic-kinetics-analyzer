"""Regenerate paper/figure.png.

(a) Global fit of the competitive-inhibition example (data/example_inhibition.csv).
(b) Coverage of nominal 95% confidence intervals for Km over 1,000 simulated experiments
per design, with intervals from the normal distribution (z = 1.96) or from Student's t
distribution with n - 2 degrees of freedom (as used by the package). Designs use n
substrate concentrations spaced logarithmically from 0.25 Km to 16 Km, with 5% constant
noise relative to Vmax.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import enzyme_kinetics as ek  # noqa: E402

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.7), dpi=300)

data = pd.read_csv(ROOT / "data" / "example_inhibition.csv")
fit = ek.fit_inhibition(data["substrate"], data["inhibitor"], data["velocity"], "competitive")
p = fit.params
grid = np.logspace(np.log10(0.4), np.log10(80), 200)
colors = plt.cm.viridis(np.linspace(0, 0.85, data["inhibitor"].nunique()))
for c, (i, g) in zip(colors, data.groupby("inhibitor")):
    m = g.groupby("substrate")["velocity"].mean()
    ax1.plot(m.index, m.values, "o", color=c, ms=4, label=f"[I] = {i:g}")
    ax1.plot(grid, ek.competitive(grid, i, p["Vmax"], p["Km"], p["Ki"]), "-", color=c, lw=1)
ax1.set_xscale("log")
ax1.set_xlabel("Substrate concentration")
ax1.set_ylabel("Initial rate")
ax1.set_title(f"(a) Competitive fit: Ki = {p['Ki']:.2f} ({fit.ci_low['Ki']:.2f} to {fit.ci_high['Ki']:.2f})",
              fontsize=9, loc="left")
ax1.legend(fontsize=7, frameon=False)

rng = np.random.default_rng(0)
ns = [4, 5, 6, 8, 10, 12, 16]
cov_z, cov_t = [], []
for n in ns:
    s = 4.0 * np.logspace(np.log10(0.25), np.log10(16), n)
    hz = ht = 0
    for _ in range(1000):
        v = ek.michaelis_menten(s, 10.0, 4.0) + rng.normal(0, 0.5, n)
        r = ek.fit_michaelis_menten(s, v)
        se = r.stderr["Km"]
        hz += abs(r.params["Km"] - 4.0) <= 1.96 * se
        ht += abs(r.params["Km"] - 4.0) <= stats.t.ppf(0.975, n - 2) * se
    cov_z.append(hz / 1000)
    cov_t.append(ht / 1000)
ax2.plot(ns, cov_z, "s--", color="#c0504d", label="Normal (z = 1.96), version 1")
ax2.plot(ns, cov_t, "o-", color="#1f6fb2", label="Student's t, version 2")
ax2.axhline(0.95, color="k", ls=":", lw=1)
ax2.set_xlabel("Substrate concentrations per experiment")
ax2.set_ylabel("Coverage of 95% interval for Km")
ax2.set_title("(b) Interval coverage in simulated experiments", fontsize=9, loc="left")
ax2.legend(fontsize=7, frameon=False, loc="lower right")
for a in (ax1, ax2):
    a.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(ROOT / "paper" / "figure.png")
print("n:", ns)
print("coverage z:", cov_z)
print("coverage t:", cov_t)
