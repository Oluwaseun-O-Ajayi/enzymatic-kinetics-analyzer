import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.dose_response import fit_ic50, logistic4
import os

# -------------------------------
# Load batch IC50 data
# -------------------------------
data = pd.read_csv("data/batch_ic50.csv")

# Create results folder if it doesn't exist
os.makedirs("results", exist_ok=True)

# Store batch results
batch_results = []

# Loop over each unique inhibitor
for inhibitor_name in data['inhibitor'].unique():
    subset = data[data['inhibitor'] == inhibitor_name]
    
    # Average replicates if multiple rows per inhibitor
    response_cols = subset.columns[1:]
    mean_response = subset[response_cols].mean(axis=0)
    std_response = subset[response_cols].std(axis=0)
    
    concentrations = np.array([float(col.replace("rep", "")) for col in response_cols])  # optional: map reps to conc
    concentrations = np.arange(1, len(response_cols)+1)  # placeholder if needed

    # Fit 4PL model
    popt = fit_ic50(concentrations, mean_response)
    bottom, top, ic50, hill = popt
    
    # Save results
    batch_results.append({
        "Inhibitor": inhibitor_name,
        "Bottom": bottom,
        "Top": top,
        "IC50": ic50,
        "Hill": hill
    })
    
    # Plot
    x_fit = np.logspace(np.log10(concentrations.min()/2), np.log10(concentrations.max()*2), 200)
    y_fit = logistic4(x_fit, *popt)
    
    plt.figure(figsize=(6, 4))
    plt.errorbar(concentrations, mean_response, yerr=std_response, fmt='o', capsize=3)
    plt.plot(x_fit, y_fit, '-', color='red')
    plt.xscale('log')
    plt.xlabel("Inhibitor concentration (µM)")
    plt.ylabel("Response (%)")
    plt.title(f"{inhibitor_name} IC50 Curve")
    plt.tight_layout()
    plt.savefig(f"figures/{inhibitor_name}_ic50.png", dpi=300)
    plt.close()

# Save batch results
results_df = pd.DataFrame(batch_results)
results_df.to_csv("results/batch_ic50_results.csv", index=False)
print("✅ Batch IC50 analysis complete. Results saved to results/batch_ic50_results.csv")
