import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.dose_response import fit_ic50, logistic4  # logistic4 = 4PL function

# -------------------------------
# Load IC50 data with replicates
# -------------------------------
data = pd.read_csv("data/example_ic50.csv")

# Compute mean and standard deviation across replicate columns
response_cols = data.columns[1:]  # all columns except inhibitor concentration
mean_response = data[response_cols].mean(axis=1)
std_response = data[response_cols].std(axis=1)

# Inhibitor concentrations
inhibitor = data['inhibitor_concentration_uM']

# -------------------------------
# Fit 4-parameter logistic (4PL) model
# -------------------------------
# popt = [bottom, top, IC50, hill_slope]
popt = fit_ic50(inhibitor, mean_response)
bottom, top, ic50, hill = popt

print("===== IC50 Fitting Results =====")
print(f"Bottom = {bottom:.2f}")
print(f"Top = {top:.2f}")
print(f"IC50 = {ic50:.3f} µM")
print(f"Hill slope = {hill:.2f}")
print("================================")

# -------------------------------
# Generate fitted curve for plotting
# -------------------------------
x_fit = np.logspace(np.log10(inhibitor.min()/2), np.log10(inhibitor.max()*2), 200)
y_fit = logistic4(x_fit, *popt)  # <-- use fitted parameters

# -------------------------------
# Plot data with error bars
# -------------------------------
plt.figure(figsize=(6, 4))
plt.errorbar(inhibitor, mean_response, yerr=std_response, fmt='o', label='Data ± SD', capsize=3)
plt.plot(x_fit, y_fit, '-', color='red', label='4PL Fit')
plt.xscale('log')
plt.xlabel("Inhibitor concentration (µM)")
plt.ylabel("Response (%)")
plt.title("IC50 Dose-Response Curve")
plt.legend()
plt.tight_layout()

# Save figure
plt.savefig("figures/ic50_with_errorbars.png", dpi=300)
plt.show()
