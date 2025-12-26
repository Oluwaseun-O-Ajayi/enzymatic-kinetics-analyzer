import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from src.dose_response import fit_ic50, logistic4  # use logistic4 for plotting

# -------------------------------
# Load batch IC50 data with replicates
# -------------------------------
data = pd.read_csv("data/batch_ic50_with_reps.csv")
drugs = data['inhibitor']

results_list = []

for idx, row in data.iterrows():
    drug = row['inhibitor']
    response_vals = row[1:].values.astype(float)
    conc_cols = row.index[1:].values.astype(str)
    
    # Extract unique concentrations from column names
    concentrations = sorted(set([float(c.split('_')[0]) for c in conc_cols]))
    
    mean_response = []
    std_response = []
    
    # Compute mean ± SD per concentration
    for conc in concentrations:
        # mask columns for this concentration
        mask = [float(c.split('_')[0]) == conc for c in conc_cols]
        rep_vals = response_vals[mask]
        mean_response.append(rep_vals.mean())
        std_response.append(rep_vals.std())
    
    mean_response = np.array(mean_response)
    std_response = np.array(std_response)
    
    # Fit 4PL IC50
    popt = fit_ic50(np.array(concentrations), mean_response)
    bottom, top, ic50, hill = popt
    results_list.append([drug, bottom, top, ic50, hill])
    
    # -------------------------------
    # Generate smooth curve for plotting using logistic4
    # -------------------------------
    x_fit = np.logspace(np.log10(min(concentrations)/2), np.log10(max(concentrations)*2), 200)
    y_fit = logistic4(x_fit, bottom, top, ic50, hill)  # directly use fitted params
    
    # Plot
    plt.figure(figsize=(6,4))
    plt.errorbar(concentrations, mean_response, yerr=std_response, fmt='o', label='Data ± SD', capsize=3)
    plt.plot(x_fit, y_fit, '-', color='red', label='4PL Fit')
    plt.xscale('log')
    plt.xlabel("Inhibitor concentration (µM)")
    plt.ylabel("Response (%)")
    plt.title(f"{drug} IC50 Dose-Response")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"figures/{drug}_ic50_replicates.png", dpi=300)
    plt.close()

# Save batch results
results_df = pd.DataFrame(results_list, columns=['Inhibitor', 'Bottom', 'Top', 'IC50', 'Hill'])
results_df.to_csv("results/batch_ic50_with_reps_results.csv", index=False)

print("✅ Batch IC50 analysis with replicates complete. Results saved to results/batch_ic50_with_reps_results.csv")
