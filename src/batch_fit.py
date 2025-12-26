import pandas as pd
from src.fit_mm import fit_michaelis_menten

def batch_fit_mm(csv_path):
    df = pd.read_csv(csv_path)
    results = []

    grouped = df.groupby(["enzyme", "condition"])

    for (enzyme, condition), group in grouped:
        vmax, km, pcov = fit_michaelis_menten(
            group["substrate_concentration_uM"],
            group["initial_velocity"]
        )

        results.append({
            "enzyme": enzyme,
            "condition": condition,
            "Vmax": vmax,
            "Km": km
        })

    return pd.DataFrame(results)
