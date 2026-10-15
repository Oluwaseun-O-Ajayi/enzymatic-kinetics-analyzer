# Example data

All files are illustrative values for demonstrating input formats and outputs. They are not experimental measurements of specific enzymes or compounds.

| File | Layout | Read with |
|---|---|---|
| `example_michaelis_menten.csv` | `substrate, velocity` | `pandas.read_csv` |
| `batch_kinetics.csv` | `enzyme, condition, substrate_concentration_uM, initial_velocity` (long format); simulated, true parameters listed in `examples/make_example_data.py` | `pandas.read_csv`, then `batch_michaelis_menten` |
| `example_ic50.csv` | first column concentration, then one column per replicate | `read_replicate_table` |
| `batch_ic50.csv` | one row per compound, one column per concentration | `read_compound_table` |
| `batch_ic50_with_reps.csv` | one row per compound, columns `concentration_RepN`; concentrations stop at 0.5, so some IC50 values lie outside the tested range | `read_compound_table` |
| `example_inhibition.csv` | `substrate, inhibitor, replicate, velocity`; simulated competitive inhibition with Vmax 12, Km 4, Ki 2.5 (see `examples/make_example_data.py`) | `pandas.read_csv`, then `compare_inhibition_models` |
