# Changelog

## 2.0.0 (2026-10)

Version 2 replaces the script collection of version 1 with an installable package, `enzyme_kinetics`. The scripts in the repository root and the `src/` folder were removed; their tasks are covered by `examples/`.

### Fixed
- `run_analysis.py` failed on the bundled data (it expected a column named `substrate_concentration_uM`).
- `run_ic50_analysis_with_reps.py` failed (it passed fitted parameters back into the fitting function as data) and averaged across concentrations rather than replicates.
- `generate_pdf_report.py` required `fpdf`, which was not listed as a dependency; reports are now built with Matplotlib.
- Confidence intervals used z = 1.96 regardless of the requested level and sample size. They now use Student's t distribution with n - p degrees of freedom; in simulated eight-point experiments the coverage of nominal 95% intervals for Km rose from about 0.89 to about 0.94.
- Two Michaelis-Menten fitters and two inconsistent IC50 models (one with the Hill slope fixed at 1, one bounded to 0 to 100) were merged into one implementation each.
- The bundled batch kinetics file had one to three points per group; it now holds simulated data with eight substrate concentrations per group and known parameters.

### Added
- IC50 fitted on the log scale with an asymmetric interval; optional fixed bottom and top.
- Replicates fitted as individual points; readers for replicate and compound tables.
- Competitive, uncompetitive and mixed inhibition models with global fitting and AICc model comparison.
- Warnings for dose-response fits that extrapolate plateaus, place the IC50 outside the tested range or have too few concentrations on one side of it.
- Tests, including a simulation check of interval coverage, and continuous integration on Linux, macOS and Windows; `CONTRIBUTING.md` and this changelog.

## 1.0.0 (2026-07)
- Initial archived release.
