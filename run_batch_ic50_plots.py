import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.dose_response import fit_ic50

# Load batch IC50 data with replicates
data = pd.read_csv("results/batch_ic50_results_with_replicates.csv")  # updated CSV

inhibitors = data['Inhibitor']
replicate_cols = data.columns[1:]  # all replicate columns

plt.figure(figsize=(8,6))

for i, row in data.iterrows():
    rep_values = row[replicate_cols].values.astype(float)
    mean_response = np.mean(rep_values)
    std_response = np.std(rep_values)

    # Inhibitor concentrations (same for all, you can store separately if needed)
    x = np.array([0.01,0.03,0.1,0.3,1,3,10])  # adjust to your experimental points
    
    # Fit 4PL curve using mean
    popt = fit_ic50(x, rep_values)
    bottom, top, ic50, hill = popt

    x_fit = np.logspace(np.log10(x.min()/2), np.log10(x.max()*2), 200)
    y_fit = fit_ic50(x_fit, rep_values, return_curve=True)

    # Plot with error bars
    plt.errorbar(x, rep_values, yerr=std_response, fmt='o', label=f"{row['Inhibitor']} ± SD", capsize=3)
    plt.plot(x_fit, y_fit, '-', label=f"{row['Inhibitor']} fit")

plt.xscale('log')
plt.xlabel("Inhibitor concentration (µM)")
plt.ylabel("Response (%)")
plt.title("Batch IC50 Curves with Replicates")
plt.legend()
plt.tight_layout()
plt.savefig("figures/batch_ic50_curves_with_errorbars.png", dpi=300)
plt.show()
