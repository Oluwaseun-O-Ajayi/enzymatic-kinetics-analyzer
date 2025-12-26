# Enzymatic Kinetics Analyzer ⚗️📊

**Python-based toolkit for comprehensive enzymatic kinetics and IC50 analysis, designed for bioanalytical chemistry, enzymology, and pharmaceutical research workflows.**

This toolkit provides robust, reproducible analysis pipelines for both single-enzyme and batch experiments, with a focus on scalability, precision, and real-world applicability. It is intended for academic research, pharmaceutical screening, and industrial biotechnology applications.


## 🚀 Key Features

* **Nonlinear Michaelis–Menten fitting** for accurate estimation of Vmax and Km  
* **4-Parameter Logistic (4PL) IC50 analysis**, including support for replicates  
* **Confidence intervals and error estimation** for meaningful statistical interpretation  
* **Batch analysis across multiple enzymes, conditions, or inhibitors** for high-throughput workflows  
* **Publication-quality plot generation** for quick visualization of dose-response curves  
* **Extensible Python framework** to integrate additional models, metrics, or automated reporting  
* **Replicates support** with automatic calculation of mean ± standard deviation for reliable experimental interpretation  


## 🎯 Motivation

This toolkit was developed to **streamline enzyme kinetics and IC50 analysis**, minimizing manual calculation errors, enabling reproducible batch processing, and producing outputs ready for publication or industry use.



## 🧩 Key Concepts

### Michaelis–Menten Kinetics

Michaelis–Menten kinetics describe the relationship between substrate concentration and enzyme-catalyzed reaction rate.

* **Vmax (Maximum Velocity):** Maximum reaction rate under saturating substrate  
* **Km (Michaelis Constant):** Substrate concentration at half-maximal velocity; indicates enzyme affinity  

**Why it matters:**  

* Quantifies enzyme efficiency and substrate preference  
* Compares enzyme variants or conditions  
* Supports mechanistic understanding and drug screening  


### IC50 and Dose-Response

**IC50** is the inhibitor concentration that reduces enzyme activity by 50%.

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



### Replicates and Statistical Confidence

* Mean ± standard deviation and confidence intervals included  
* Facilitates reproducible comparisons across experimental conditions or compounds  



## 📁 Folder Structure



enzymatic-kinetics-analyzer/
├── data/                        # Example CSV datasets
├── figures/                     # Generated plots for visualization
├── results/                     # Output tables and reports
├── src/                          # Core Python functions
├── run_analysis.py               # Single enzyme kinetic analysis
├── run_batch_analysis.py         # Batch Michaelis–Menten analysis
├── run_ic50_analysis.py          # Single IC50 analysis
├── run_batch_ic50_analysis.py    # Batch IC50 analysis
├── run_batch_ic50_with_reps.py   # Batch IC50 with replicates
├── generate_pdf_report.py        # PDF report generator
├── requirements.txt              # Python dependencies
└── README.md



## ⚡ Quick-start Snippet

Run complete batch analysis and generate PDF report:


python run_batch_analysis.py && python generate_pdf_report.py




## 🧪 Workflow Overview


          ┌───────────────┐
          │ Prepare CSV   │
          │ - Michaelis–  │
          │   Menten or   │
          │   IC50 data   │
          │ - Include     │
          │   replicates  │
          └───────┬───────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Input Validation    │
       │ - Check columns     │
       │ - Ensure numeric    │
       │   data              │
       │ - Missing values    │
       │   handled           │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Optional Preprocess │
       │ - Normalize data    │
       │ - Correct units     │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Choose Analysis     │
       │ - Single enzyme MM  │
       │ - Batch MM          │
       │ - Single IC50       │
       │ - Batch IC50 w/ reps│
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Run Script           │
       │ e.g., run_analysis.py│
       │      run_batch_      │
       │      analysis.py     │
       │      run_ic50_       │
       │      analysis.py     │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Computation & Fitting│
       │ - Michaelis–Menten   │
       │   Vmax & Km          │
       │ - IC50 (4PL)         │
       │ - Mean ± SD for reps │
       │ - Confidence intervals│
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Review Plots & QC   │
       │ - Check curve fits  │
       │ - Identify outliers │
       │ - Validate results  │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Save Outputs        │
       │ - CSV tables        │
       │ - Figures (plots)   │
       │ - Batch summary     │
       │   tables            │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Generate PDF Report │
       │ - Combines plots &  │
       │   tables            │
       │ - Unicode support   │
       │ - Publication-ready │
       └─────────┬───────────┘
                  │
                  ▼
       ┌─────────────────────┐
       │ Downstream Analysis  │
       │ - Compare enzymes &  │
       │   inhibitors         │
       │ - Statistical insights│
       │ - Data-driven        │
       │   conclusions        │
       └─────────────────────┘


## 📝 Example CSV Formats

**Michaelis–Menten (single enzyme):**


substrate_concentration_uM,initial_velocity
0.5,0.08
1.0,0.15
2.0,0.29
5.0,0.63
10.0,0.90


**IC50 with replicates (batch):**


inhibitor,0.01_Rep1,0.01_Rep2,0.01_Rep3,0.05_Rep1,0.05_Rep2,0.05_Rep3,0.1_Rep1,0.1_Rep2,0.1_Rep3
DrugA,95,96,94,91,92,90,84,85,83
DrugB,90,88,91,89,90,87,85,84,86
DrugC,85,87,86,83,82,84,80,81,79


## 💻 Installation


pip install -r requirements.txt


## Usage

### Single Analysis


python run_analysis.py


* Michaelis–Menten fit for a single enzyme dataset
* Outputs Vmax, Km, and plots in `figures/`

### Batch Analysis

python run_batch_analysis.py

* Processes multiple enzymes or conditions
* Aggregates results in `results/batch_kinetics_results.csv`

### IC50 Analysis

python run_ic50_analysis.py
python run_batch_ic50_analysis.py
python run_batch_ic50_with_reps.py


* Fits 4PL dose-response curves
* Supports replicates and batch analysis
* Saves plots and aggregated results in `figures/` and `results/`


## 📊 Output

* CSV tables of kinetic and IC50 parameters
* Error bars, replicate statistics, and confidence intervals included
* Publication-quality **plots**: Michaelis–Menten curves and IC50 dose-response curves


## 📂 Output Files

* `results/batch_kinetics_results.csv`
* `results/batch_ic50_with_reps_results.csv`
* `results/enzymatic_kinetics_report.pdf`


## ⚠️ Note About Fonts

Ensure `DejaVuSans.ttf` is in `fonts/` and not read-only for proper Unicode support in PDF reports.


## 🌐 Professional Context

This toolkit is designed to **bridge research and industry needs**:

* Automated, reproducible workflows for **high-throughput enzyme screening**
* Supports **drug discovery pipelines** with rigorous IC50 and kinetic analysis
* Facilitates **data-driven decision-making** with clear outputs, plots, and replicate management
* Extensible for integration into **automated lab systems** or proprietary analytics pipelines

Positions the user as a **valuable contributor to bioanalytical, pharmaceutical, or biotech teams**, combining technical expertise, analytical rigor, and software-driven efficiency


## ✍️ Author

**Oluwaseun O. Ajayi**
PhD Researcher, Chemistry
Bioanalytical & Structural Chemistry | Enzymology | Lab Automation
