from src.batch_fit import batch_fit_mm

results = batch_fit_mm("data/batch_kinetics.csv")
results.to_csv("results/batch_kinetics_results.csv", index=False)

print(results)
