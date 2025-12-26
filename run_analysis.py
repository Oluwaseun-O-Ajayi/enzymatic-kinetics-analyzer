import pandas as pd
import matplotlib.pyplot as plt

from src.michaelis_menten import fit_michaelis_menten, michaelis_menten
from src.fitting_utils import parameter_confidence_intervals

data = pd.read_csv("data/example_michaelis_menten.csv")
  

data = pd.read_csv("data/example_michaelis_menten.csv")

substrate = data["substrate_concentration_uM"].astype(float)
velocity = data["initial_velocity"].astype(float)



popt, pcov = fit_michaelis_menten(substrate, velocity)
vmax, km = popt

intervals = parameter_confidence_intervals(popt, pcov)

print(f"Vmax = {vmax:.3f} µM/min")
print(f"Km = {km:.3f} µM")
print("Confidence intervals:", intervals)

import numpy as np

substrate_fit = np.array(sorted(substrate))

velocity_fit = michaelis_menten(substrate_fit, vmax, km)

plt.scatter(substrate, velocity)
plt.plot(substrate_fit, velocity_fit)
plt.xlabel("Substrate concentration (µM)")
plt.ylabel("Initial velocity (µM/min)")
plt.tight_layout()
plt.savefig("figures/michaelis_menten_plot.png", dpi=300)
plt.show()

