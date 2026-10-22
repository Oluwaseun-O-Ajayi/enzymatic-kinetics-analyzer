"""Tests for Enzymatic Kinetics Analyzer: parameter recovery, interval coverage, model
selection, input readers and example scripts."""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import enzyme_kinetics as ek

ROOT = Path(__file__).resolve().parents[1]
S = np.array([0.5, 1, 2, 4, 8, 16, 32, 64.0])


def test_michaelis_menten_recovers_exact_parameters():
    fit = ek.fit_michaelis_menten(S, ek.michaelis_menten(S, 12.0, 4.0))
    assert fit.params["Vmax"] == pytest.approx(12.0, rel=1e-6)
    assert fit.params["Km"] == pytest.approx(4.0, rel=1e-6)
    assert fit.r_squared == pytest.approx(1.0)


def test_confidence_interval_uses_t_distribution():
    rng = np.random.default_rng(0)
    v = ek.michaelis_menten(S, 12.0, 4.0) + rng.normal(0, 0.2, S.size)
    fit = ek.fit_michaelis_menten(S, v)
    from scipy import stats
    t = stats.t.ppf(0.975, S.size - 2)
    assert fit.ci_high["Km"] - fit.params["Km"] == pytest.approx(t * fit.stderr["Km"])
    narrow = ek.fit_michaelis_menten(S, v, confidence=0.80)
    assert narrow.ci_high["Km"] - narrow.ci_low["Km"] < fit.ci_high["Km"] - fit.ci_low["Km"]


def test_km_interval_coverage_is_close_to_nominal():
    """Over repeated simulated experiments, the 95% interval for Km should contain the
    true value close to 95% of the time when the design brackets Km."""
    rng = np.random.default_rng(1)
    hits = 0
    n_sim = 300
    for _ in range(n_sim):
        v = ek.michaelis_menten(S, 12.0, 4.0) + rng.normal(0, 0.3, S.size)
        fit = ek.fit_michaelis_menten(S, v)
        hits += fit.ci_low["Km"] <= 4.0 <= fit.ci_high["Km"]
    assert 0.90 <= hits / n_sim <= 0.99


def test_replicates_are_all_fitted():
    s = np.repeat(S, 3)
    v = ek.michaelis_menten(s, 12.0, 4.0) + np.random.default_rng(2).normal(0, 0.2, s.size)
    fit = ek.fit_michaelis_menten(s, v)
    assert fit.n == 24 and fit.dof == 22


def test_dose_response_recovers_ic50_and_hill():
    x = np.logspace(-3, 2, 11)
    y = ek.four_parameter_logistic(x, 5.0, 95.0, np.log10(0.8), 1.2)
    fit = ek.fit_dose_response(x, y)
    assert fit.derived["IC50"]["value"] == pytest.approx(0.8, rel=1e-4)
    assert fit.params["hill"] == pytest.approx(1.2, rel=1e-4)


def test_ic50_interval_is_asymmetric_and_contains_estimate():
    x = np.repeat(np.logspace(-2, 2, 9), 3)
    y = ek.four_parameter_logistic(x, 0.0, 100.0, 0.0, 1.0) + np.random.default_rng(3).normal(0, 3, x.size)
    ic = ek.fit_dose_response(x, y).derived["IC50"]
    assert ic["ci_low"] < ic["value"] < ic["ci_high"]
    assert ic["value"] / ic["ci_low"] == pytest.approx(ic["ci_high"] / ic["value"], rel=1e-6)


def test_fixed_bottom_and_top():
    x = np.logspace(-2, 2, 9)
    y = ek.four_parameter_logistic(x, 0.0, 100.0, np.log10(2.0), 1.0)
    fit = ek.fit_dose_response(x, y, fix_bottom=0.0, fix_top=100.0)
    assert fit.params["bottom"] == 0.0 and fit.params["top"] == 100.0
    assert fit.derived["IC50"]["value"] == pytest.approx(2.0, rel=1e-4)
    assert fit.dof == x.size - 2


def test_dose_response_rejects_zero_concentration():
    with pytest.raises(ValueError):
        ek.fit_dose_response([0, 1, 2, 3, 4], [100, 80, 50, 20, 5])


@pytest.mark.parametrize("truth", ["competitive", "uncompetitive"])
def test_model_comparison_identifies_the_simulated_mechanism(truth):
    rng = np.random.default_rng(4)
    s = np.tile(S, 4)
    i = np.repeat([0.0, 2.0, 5.0, 10.0], S.size)
    v = getattr(ek, truth)(s, i, 12.0, 4.0, 2.5) * (1 + rng.normal(0, 0.02, s.size))
    ranking = ek.compare_inhibition_models(s, i, v, models_to_fit=("competitive", "uncompetitive"))
    assert ranking.loc[0, "model"] == truth
    assert ranking["akaike_weight"].sum() == pytest.approx(1.0)
    fit = ek.fit_inhibition(s, i, v, truth)
    assert fit.params["Ki"] == pytest.approx(2.5, rel=0.1)


def test_mixed_model_reduces_to_noncompetitive_when_alpha_is_one():
    s, i = np.array([1.0, 4.0]), np.array([2.0, 2.0])
    v_mixed = ek.mixed(s, i, 12.0, 4.0, 2.0, 1.0)
    v_noncomp = 12.0 * s / ((4.0 + s) * (1 + i / 2.0))
    np.testing.assert_allclose(v_mixed, v_noncomp)


def test_readers_produce_long_format():
    reps = ek.read_replicate_table(ROOT / "data" / "example_ic50.csv")
    assert list(reps.columns) == ["concentration", "replicate", "response"]
    raw = pd.read_csv(ROOT / "data" / "example_ic50.csv")
    assert len(reps) == raw.shape[0] * (raw.shape[1] - 1)
    long = ek.read_compound_table(ROOT / "data" / "batch_ic50_with_reps.csv")
    assert set(long["replicate"]) == {"Rep1", "Rep2", "Rep3"}
    assert long.groupby("compound").size().eq(12).all()


def test_batch_dose_response_flags_ic50_outside_range():
    long = ek.read_compound_table(ROOT / "data" / "batch_ic50.csv")
    table = ek.batch_dose_response(long)
    assert len(table) == 3 and (table["note"] == "").all()


@pytest.mark.parametrize("script", ["michaelis_menten.py", "dose_response.py", "inhibition.py", "report.py"])
def test_examples_run(script, tmp_path):
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run([sys.executable, str(ROOT / "examples" / script)], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]


@pytest.mark.parametrize("script", sorted(p.name for p in (ROOT / "examples").glob("case_study_*.py")))
def test_case_studies_run(script, tmp_path):
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.run([sys.executable, str(ROOT / "examples" / script)], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-2000:]
