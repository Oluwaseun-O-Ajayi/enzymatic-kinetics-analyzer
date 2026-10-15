"""Readers that convert common spreadsheet layouts to long format.

Long format has one row per measurement, which keeps replicates as separate points so
that every measurement enters the fit.
"""

import re
from pathlib import Path
from typing import Union

import pandas as pd

PathLike = Union[str, Path]


def read_replicate_table(path: PathLike) -> pd.DataFrame:
    """Read a table whose first column is concentration and whose other columns are
    replicate responses (for example `conc, rep1, rep2, rep3`).

    Returns columns `concentration`, `replicate` and `response`."""
    df = pd.read_csv(path)
    conc_col = df.columns[0]
    long = df.melt(id_vars=conc_col, var_name="replicate", value_name="response")
    long = long.rename(columns={conc_col: "concentration"}).dropna(subset=["response"])
    return long[["concentration", "replicate", "response"]].reset_index(drop=True)


_HEADER = re.compile(r"^\s*([0-9.eE+-]+)(?:_(.+))?\s*$")


def read_compound_table(path: PathLike, id_column: str = None) -> pd.DataFrame:
    """Read a table with one row per compound and one column per concentration.

    Column headers are concentrations, optionally followed by a replicate label after an
    underscore (`0.1`, or `0.1_Rep1`, `0.1_Rep2`, ...). The first column (or `id_column`)
    identifies the compound.

    Returns columns `compound`, `concentration`, `replicate` and `response`."""
    df = pd.read_csv(path)
    id_column = id_column or df.columns[0]
    rows = []
    for col in df.columns:
        if col == id_column:
            continue
        m = _HEADER.match(str(col))
        if not m:
            raise ValueError(f"Cannot read a concentration from column header '{col}'")
        conc, rep = float(m.group(1)), m.group(2) or "1"
        for compound, value in zip(df[id_column], df[col]):
            if pd.notna(value):
                rows.append({"compound": compound, "concentration": conc,
                             "replicate": rep, "response": float(value)})
    return pd.DataFrame(rows, columns=["compound", "concentration", "replicate", "response"])
