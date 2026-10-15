# Enzymatic Kinetics Analyzer

**Michaelis-Menten, dose-response and reversible-inhibition fitting in Python, with t-based confidence intervals, replicate-aware batch processing and model comparison.**

[![tests](https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer/actions/workflows/tests.yml/badge.svg)](https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer/actions/workflows/tests.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21652522.svg)](https://doi.org/10.5281/zenodo.21652522)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

The package fits the standard models of steady-state enzyme kinetics and inhibitor potency to tabular data and reports every parameter with its standard error and confidence interval. Replicates are kept as separate points, so they enter the fit and the degrees of freedom; spreadsheet layouts common in the laboratory are converted to long format by two readers; and results for many enzymes, conditions or compounds are collected in one table.

## What it does

| Task | Function | Notes |
|---|---|---|
| Michaelis-Menten | `fit_michaelis_menten` | Optional weights from measurement standard deviations |
| Dose-response (IC50) | `fit_dose_response` | Four-parameter logistic fitted on log10(IC50), so the IC50 interval is asymmetric; bottom and top can be fixed |
| Reversible inhibition | `fit_inhibition` | Competitive, uncompetitive and mixed (noncompetitive when alpha = 1), fitted globally to all substrate and inhibitor concentrations |
| Model comparison | `compare_inhibition_models` | Ranks models by small-sample AIC (AICc) and reports Akaike weights |
| Batch fitting | `batch_michaelis_menten`, `batch_dose_response` | One row per group; dose-response fits carry plain-language warnings when the design does not support the estimate |
| Input | `read_replicate_table`, `read_compound_table` | Concentration-by-replicate tables and compound-by-concentration tables (headers such as `0.1` or `0.1_Rep2`) |
| Output | `plotting`, `report.write_pdf_report` | Figures and a multi-page PDF built with Matplotlib only |

Confidence intervals use Student's t distribution with n - p degrees of freedom. They are asymptotic intervals: when the data constrain a parameter poorly (for example, Km when no substrate concentration exceeds Km), they can be too narrow, and a better experimental design is the remedy.

## Installation

Requires Python 3.9 or later.

```bash
git clone https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer.git
cd enzymatic-kinetics-analyzer
pip install .            # or: pip install ".[test]" to include the test dependencies
```

## Quick start

```python
import numpy as np
import enzyme_kinetics as ek

s = np.array([0.5, 1, 2, 4, 8, 16, 32, 64])
v = np.array([1.36, 2.42, 4.11, 5.98, 7.79, 9.62, 10.75, 11.17])
fit = ek.fit_michaelis_menten(s, v)
print(fit.table())          # Vmax and Km with standard errors and 95% intervals

reps = ek.read_replicate_table("data/example_ic50.csv")
ic = ek.fit_dose_response(reps["concentration"], reps["response"])
print(ic.derived["IC50"])  # value and asymmetric 95% interval
```

Worked examples: `examples/michaelis_menten.py`, `examples/dose_response.py`, `examples/inhibition.py` and `examples/report.py`. The bundled data are illustrative or simulated with known parameters (see `data/README.md`); they are not measurements of specific enzymes or compounds.

## Testing

```bash
pip install ".[test]"
python -m pytest
```

The tests check exact recovery of parameters from noise-free data, that intervals use the t distribution, that 95% intervals for Km cover the true value close to 95% of the time in repeated simulated experiments, that replicates enter the fit, that the IC50 interval is symmetric on the log scale, fixed plateaus, the identification of the simulated inhibition mechanism, the input readers, and that all examples run. Tests run on Linux, macOS and Windows for every push.

## Contributing and support

Bug reports, questions and pull requests are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md). Please use [GitHub Issues](https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer/issues); the author aims to respond within two weeks.

## Citation

See `CITATION.cff`.

## License

MIT. See [LICENSE](LICENSE).

## Support

Questions, bug reports and feature requests are welcome as GitHub issues. The author maintains the package and reviews issues and pull requests; see CONTRIBUTING.md for how to contribute.
