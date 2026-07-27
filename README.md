# Enzymatic Kinetics Analyzer

**Open-source Python toolkit for reproducible enzymatic kinetics and dose–response analysis**

**The framework provides Michaelis–Menten parameter estimation, 4-parameter logistic (4PL) IC50 modeling, confidence interval calculation, batch processing workflows, publication-quality visualization, and report generation for enzymology, bioanalytical chemistry, and pharmaceutical research applications**

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Enzymatic Kinetics Analyzer is a reproducible computational workflow for analyzing enzyme kinetics and inhibitor potency data

The software integrates:

- Michaelis–Menten kinetic modeling
- IC50 and dose-response analysis
- Confidence interval estimation
- Replicate-aware statistical analysis
- Batch processing workflows
- Publication-quality figure generation
- Automated PDF reporting

The framework is intended for educational, academic, industrial, and pharmaceutical research environments where transparent and reproducible enzyme characterization workflows are required

## Scope and Intended Use

This repository provides computational workflows for enzyme kinetic characterization and inhibitor potency analysis

The software performs:

- Model fitting
- Statistical analysis
- Visualization
- Reporting

The repository does not replace experimental validation or laboratory quality systems

Users are responsible for validating fitted parameters, assay quality, and experimental assumptions before applying results to regulated or decision-critical workflows

## Key Features

* **Nonlinear Michaelis–Menten fitting** for accurate estimation of Vmax and Km  
* **4-Parameter Logistic (4PL) IC50 analysis**, including support for replicates  
* **Confidence intervals and error estimation** for meaningful statistical interpretation  
* **Batch analysis across multiple enzymes, conditions, or inhibitors** for high-throughput workflows  
* **Publication-quality plot generation** for quick visualization of dose-response curves  
* **Extensible Python framework** to integrate additional models, metrics, or automated reporting  
* **Replicates support** with automatic calculation of mean ± standard deviation for reliable experimental interpretation  

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer.git
cd enzymatic-kinetics-analyzer

# Install dependencies
pip install -r requirements.txt
```

### Quick Analysis

```bash
# Single enzyme analysis
python run_analysis.py

# Batch analysis
python run_batch_analysis.py

# IC50 analysis
python run_ic50_analysis.py

# Generate PDF report
python generate_pdf_report.py
```

## Requirements

```
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.4.0
scipy>=1.7.0
reportlab>=3.6.0
```

## Key Concepts

### Michaelis–Menten Kinetics

Michaelis–Menten kinetics describe the relationship between substrate concentration and enzyme-catalyzed reaction rate

* **Vmax (Maximum Velocity):** Maximum reaction rate under saturating substrate  
* **Km (Michaelis Constant):** Substrate concentration at half-maximal velocity; indicates enzyme affinity  

**Why it matters:**  
* Quantifies enzyme efficiency and substrate preference  
* Compares enzyme variants or conditions  
* Supports mechanistic understanding and drug screening  

**Equation:**
```
v = (Vmax × [S]) / (Km + [S])
```

### IC50 and Dose-Response

**IC50** is the inhibitor concentration that reduces enzyme activity by 50%

* **Low IC50** → potent inhibitor  
* **High IC50** → weak inhibitor  

**Analysis approach:**  
* 4-parameter logistic (4PL) curve fitting  
* Incorporates replicates and error statistics  
* Outputs confidence intervals and publication-quality plots  

**Applications:**  
* Drug candidate screening  
* Comparative potency analysis  
* Enzyme modulation studies  

**4PL Equation:**
```
y = Bottom + (Top - Bottom) / (1 + (x/IC50)^HillSlope)
```

### Replicates and Statistical Confidence

* Mean ± standard deviation and confidence intervals included  
* Facilitates reproducible comparisons across experimental conditions or compounds  

## Project Structure

```
enzymatic-kinetics-analyzer/
├── data/                        # Example CSV datasets
├── figures/                     # Generated plots for visualization
├── results/                     # Output tables and reports
├── src/                         # Core Python functions
├── run_analysis.py              # Single enzyme kinetic analysis
├── run_batch_analysis.py        # Batch Michaelis–Menten analysis
├── run_ic50_analysis.py         # Single IC50 analysis
├── run_batch_ic50_analysis.py   # Batch IC50 analysis
├── run_batch_ic50_with_reps.py  # Batch IC50 with replicates
├── generate_pdf_report.py       # PDF report generator
├── requirements.txt             # Python dependencies
├── LICENSE                      # MIT License
├── .gitignore                   # Git ignore rules
└── README.md                    # This file
```

## Data Format

### Michaelis–Menten (single enzyme)

```csv
substrate_concentration_uM,initial_velocity
0.5,0.08
1.0,0.15
2.0,0.29
5.0,0.63
10.0,0.90
```

### IC50 with replicates (batch)

```csv
inhibitor,0.01_Rep1,0.01_Rep2,0.01_Rep3,0.05_Rep1,0.05_Rep2,0.05_Rep3,0.1_Rep1,0.1_Rep2,0.1_Rep3
DrugA,95,96,94,91,92,90,84,85,83
DrugB,90,88,91,89,90,87,85,84,86
DrugC,85,87,86,83,82,84,80,81,79
```

## Usage Examples

### Single Enzyme Analysis

```python
python run_analysis.py
```

**Output:**
* Michaelis–Menten fit for a single enzyme dataset
* Vmax and Km values with confidence intervals
* Publication-quality plots in `figures/`

### Batch Analysis

```python
python run_batch_analysis.py
```

**Output:**
* Processes multiple enzymes or conditions
* Aggregates results in `results/batch_kinetics_results.csv`
* Comparative plots for all enzymes

### IC50 Analysis

```python
# Single IC50
python run_ic50_analysis.py

# Batch IC50
python run_batch_ic50_analysis.py

# Batch IC50 with replicates
python run_batch_ic50_with_reps.py
```

**Output:**
* Fits 4PL dose-response curves
* Supports replicates and batch analysis
* IC50 values with confidence intervals
* Plots saved in `figures/` and results in `results/`

### Generate PDF Report

```python
python generate_pdf_report.py
```

**Output:**
* Comprehensive PDF report combining all analyses
* Publication-ready with plots and tables
* Unicode support for special characters

## Output Files

### Results Tables
* `results/batch_kinetics_results.csv` - Vmax, Km for all enzymes
* `results/batch_ic50_with_reps_results.csv` - IC50 values with statistics
* `results/enzymatic_kinetics_report.pdf` - Complete analysis report

### Plots
* `figures/michaelis_menten_*.png` - Kinetic curves
* `figures/ic50_*.png` - Dose-response curves
* `figures/batch_comparison.png` - Side-by-side comparisons

All outputs include:
* Error bars and confidence intervals
* Replicate statistics (mean ± SD)
* R² values for goodness of fit
* Publication-quality formatting

## Workflow Diagram

```
┌───────────────┐
│ Prepare CSV   │
│ Data          │
└───────┬───────┘
        │
        ▼
┌─────────────────────┐
│ Input Validation    │
│ - Check columns     │
│ - Numeric data      │
│ - Handle missing    │
└─────────┬───────────┘
        │
        ▼
┌─────────────────────┐
│ Choose Analysis     │
│ - Michaelis–Menten  │
│ - IC50              │
│ - Batch processing  │
└─────────┬───────────┘
        │
        ▼
┌─────────────────────┐
│ Run Analysis        │
│ - Curve fitting     │
│ - Statistics        │
│ - Error estimation  │
└─────────┬───────────┘
        │
        ▼
┌─────────────────────┐
│ Review Results      │
│ - Check plots       │
│ - Validate fits     │
│ - Identify outliers │
└─────────┬───────────┘
        │
        ▼
┌─────────────────────┐
│ Generate Report     │
│ - PDF output        │
│ - CSV tables        │
│ - Publication plots │
└─────────────────────┘
```

## Research Applications

### Drug Discovery
* **Inhibitor screening** - Rank drug candidates by IC50
* **Lead optimization** - Compare compound potency
* **Target validation** - Confirm enzyme inhibition

### Enzyme Engineering
* **Mutant comparison** - Evaluate Km and Vmax changes
* **Substrate preference** - Test multiple substrates
* **Optimization** - Improve catalytic efficiency

### Bioprocess Development
* **Enzyme selection** - Choose best catalyst
* **Process optimization** - Maximize reaction rates
* **Quality control** - Ensure batch consistency

### Academic Research
* **Mechanistic studies** - Understand enzyme function
* **Structure-function** - Relate activity to structure
* **Publication** - Generate figures and data tables

## Troubleshooting

### Fonts Issue (PDF Reports)
Ensure `DejaVuSans.ttf` is in `fonts/` directory and not read-only for proper Unicode support in PDF reports

```bash
# Check if font exists
ls fonts/DejaVuSans.ttf

# If missing, download from:
# https://dejavu-fonts.github.io/
```

### Poor Curve Fits
* **Check data quality** - Remove outliers
* **Verify concentration range** - Ensure proper range for Km or IC50
* **Increase data points** - More points = better fits
* **Check initial guesses** - Adjust parameters in code if needed

### Import Errors
```bash
# Reinstall dependencies
pip install -r requirements.txt --upgrade
```

## Customization

### Modify Fitting Parameters

Edit the analysis scripts to adjust:
* Initial parameter guesses
* Fitting bounds
* Convergence criteria
* Plot styles and colors

### Add Custom Metrics

Extend the analysis by adding:
* Kcat calculations
* Catalytic efficiency (Kcat/Km)
* Custom error models
* Additional statistical tests

## Reproducible Research

The software emphasizes transparent workflows, traceable outputs, and reproducible analysis practices

Generated outputs include:

- Figures
- CSV result tables
- Statistical summaries
- PDF reports

These outputs help researchers document and communicate kinetic and pharmacological analyses

## Contributing

Contributions welcome! Areas for enhancement:
* Additional kinetic models (e.g., competitive/non-competitive inhibition)
* Support for more file formats (Excel, JSON)
* Interactive plotting with Plotly
* Web interface for analysis
* Integration with LIMS systems

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Developed for enzymology and pharmaceutical research
* Thanks to the scientific Python community (NumPy, SciPy, Matplotlib)
* Inspired by industry-standard enzyme kinetics software

## Contact

**Oluwaseun O. Ajayi**  
PhD Researcher, Chemistry  
University of Georgia  
Bioanalytical & Structural Chemistry | Enzymology | Lab Automation

- **GitHub**: [@Oluwaseun-O-Ajayi](https://github.com/Oluwaseun-O-Ajayi)
- **Academic Email**: oluwaseun.ajayi@uga.edu
- **Personal Email**: seunolanikeajayi@gmail.com

## Publication

Associated manuscript:

Ajayi OO.

Enzymatic Kinetics Analyzer: A Reproducible Python Workflow for Michaelis–Menten and IC50 Analysis

Manuscript in preparation

## Citation

If you use this toolkit in your research:

```bibtex
@software{ajayi2026eka,
 ...
 doi={actual DOI}
}
@software{enzymatic_kinetics_analyzer,
  author = {Oluwaseun O. Ajayi},
  title = {Enzymatic Kinetics Analyzer},
  year = {2024},
  url = {https://github.com/Oluwaseun-O-Ajayi/enzymatic-kinetics-analyzer}
}
```
## Validation and Interpretation

Model quality should always be evaluated alongside experimental design, assay integrity, and goodness-of-fit metrics

Recommended practices include:

- Inspect fitted curves visually
- Review confidence intervals
- Verify concentration ranges
- Compare replicate consistency
- Confirm biological plausibility

Software outputs should support, not replace, scientific judgment

## Scientific Software Contribution

The primary contribution of this repository is a reusable software workflow for enzymatic data analysis

Implemented capabilities include:

- Nonlinear Michaelis–Menten fitting
- 4-parameter logistic dose–response fitting
- Confidence interval estimation
- Replicate-aware statistics
- Batch processing
- Automated figure generation
- Automated PDF report creation

The framework is designed to improve reproducibility and accessibility of enzymatic data analysis

---

**Advancing reproducible enzymatic data analysis through open scientific software**
