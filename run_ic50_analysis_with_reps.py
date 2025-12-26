import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.dose_response import fit_ic50

# -------------------------------
# Load IC50 data with replicates
# -------------------------------
data = pd.read_csv("data/example_ic50_with_reps.csv")

# Extract inhibitor names
drugs = data['inhibitor']

# Compute mean and standard deviation across replicates
response_cols = data.columns[1:]  # All columns except 'inhibitor'
mean_response = data[response_cols].mean(axis=1)
std_response = data[response_cols].std(axis=1)

# Example: assuming all columns represent the same concentrations
# Create an array of concentrations by extracting numbers from column names
concentrations = sorted(set(float(col.split('_')[0]) for col in response_cols))

# -------------------------------
# Fit IC50
# -------------------------------
popt = fit_ic50(np.array(concentrations), mean_response)

bottom, top, ic50, hill = popt
print("===== IC50 Fitting Results =====")
print(f"Bottom = {bottom:.2f}")
print(f"Top = {top:.2f}")
print(f"IC50 = {ic50:.3f} µM")
print(f"Hill slope = {hill:.2f}")
print("================================")

# -------------------------------
# Generate fitted curve
# -------------------------------
x_fit = np.logspace(np.log10(min(concentrations)/2), np.log10(max(concentrations)*2), 200)
y_fit = fit_ic50(x_fit, popt, return_curve=True)

# -------------------------------
# Plot data with error bars
# -------------------------------
plt.figure(figsize=(6, 4))
plt.errorbar(concentrations, mean_response, yerr=std_response, fmt='o', label='Data ± SD', capsize=3)
plt.plot(x_fit, y_fit, '-', color='red', label='4PL Fit')
plt.xscale('log')
plt.xlabel("Inhibitor concentration (µM)")
plt.ylabel("Response (%)")
plt.title("IC50 Dose-Response Curve with Replicates")
plt.legend()
plt.tight_layout()
plt.savefig("figures/ic50_with_replicates.png", dpi=300)
plt.show()
