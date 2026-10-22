# Case studies

Four worked examples that check the package against known answers. Each runs on
data bundled with the repository or on simulations with a fixed seed, writes its
results to `examples/case_study_results/`, and can be reproduced with one command:

```bash
python examples/case_study_1_parameter_recovery.py
python examples/case_study_2_inhibition_mechanism.py
python examples/case_study_3_interval_coverage.py
python examples/case_study_4_dose_response_warnings.py
```

The numbers below were produced with version 2.0.1. The simulated data sets and
their true parameters are defined in `examples/make_example_data.py`.

| Case | Application | Input |
|---|---|---|
| 1 | Recovering Vmax and Km, with honest intervals | `data/batch_kinetics.csv` (4 enzyme/condition pairs, 8 concentrations) |
| 2 | Identifying an inhibition mechanism | `data/example_inhibition.csv` (simulated competitive inhibition) |
| 3 | Checking interval coverage at small sample sizes | 2,000 simulated experiments |
| 4 | Recognizing dose-response fits that need a second look | `data/batch_ic50_with_reps.csv` (3 compounds, triplicates) |

## Case 1: Recovering known Michaelis-Menten parameters

**Question.** When the true Vmax and Km are known, how close are the fitted values,
and do the 95% confidence intervals contain the truth?

**Result.** All eight fitted parameters are within 10% of their true values. The
true value lies inside the 95% interval for 7 of 8 parameters. The exception is
Vmax for EnzB/WT (true 0.80, fitted 0.864, interval 0.831 to 0.897). The fits also
separate the wild type from the mutant: Km 10.9 (9.7 to 12.1) against 25.6 (21.9
to 29.3), with non-overlapping intervals.

**What it shows.** A 95% interval is expected to miss about 1 time in 20, so one
miss in eight is consistent with correct intervals. Reporting the miss, rather
than hiding it, is how users learn what an interval promises.

## Case 2: Identifying an inhibition mechanism

**Question.** Given rates at five inhibitor and eight substrate concentrations,
does model comparison choose the mechanism that generated the data?

**Result.** The competitive model has the lowest AICc and an Akaike weight of
0.82. The mixed model follows (delta AICc 3.1, weight 0.18) and the uncompetitive
model is ruled out (delta AICc 234). The competitive fit returns Vmax 12.00, Km
4.04 and Ki 2.59 (95% interval 2.40 to 2.79), against true values of 12, 4 and 2.5.

**What it shows.** The mixed model contains the competitive model as a limiting
case, so it fits almost as well. AICc penalizes its extra parameter, and the
weights express how strongly the data prefer one mechanism over another instead
of declaring a single winner.

## Case 3: Interval coverage at small sample sizes

**Question.** With six or eight substrate concentrations, how often does a nominal
95% interval for Km contain the true value?

**Result.** Out of 1,000 simulated experiments for each design:

| Points | Student's t (package) | Normal approximation |
|---|---|---|
| 6 | 0.908 | 0.808 |
| 8 | 0.943 | 0.884 |

**What it shows.** Intervals built on the normal distribution are too narrow when
only six to eight points are fitted. The t distribution brings coverage close to
the nominal value at eight points. At six points both fall short, because the fit
is unweighted while the simulated noise grows with the rate; adding concentrations
or weighting the fit is the remedy.

## Case 4: Dose-response fits that need a second look

**Question.** For three compounds measured in triplicate, which IC50 estimates
can be reported as they stand?

**Result.**

| Compound | IC50 (95% CI) | Hill | Note from the package |
|---|---|---|---|
| DrugA | 0.163 (0.102 to 0.262) | 1.73 | fewer than two concentrations on one side of the IC50 |
| DrugB | 0.088 (0.083 to 0.093) | 4.11 | steep Hill slope; check the concentration spacing |
| DrugC | 0.071 (0.062 to 0.081) | 0.90 | plateau extrapolated beyond the data; consider fixing bottom or top |

**What it shows.** Each fit converges and reports a tight-looking interval, yet
each has a design problem that the interval alone does not reveal. The warnings
point to the specific fix: more concentrations around the IC50, finer spacing, or
a fixed plateau.
