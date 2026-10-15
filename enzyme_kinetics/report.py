"""Multi-page PDF report built with Matplotlib (no extra dependencies)."""

from pathlib import Path
from typing import Iterable, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def write_pdf_report(path, title: str, sections: Iterable[Tuple[str, pd.DataFrame, "plt.Figure"]]):
    """Write a PDF with a title page and, for each section, a results table and an
    optional figure. `sections` is an iterable of (heading, table, figure_or_None)."""
    path = Path(path)
    with PdfPages(path) as pdf:
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(0.5, 0.6, title, ha="center", fontsize=18)
        pdf.savefig(fig)
        plt.close(fig)
        for heading, table, figure in sections:
            fig, ax = plt.subplots(figsize=(8.5, 11))
            ax.axis("off")
            ax.set_title(heading, fontsize=14, loc="left")
            shown = table.copy()
            for c in shown.select_dtypes("number").columns:
                shown[c] = shown[c].map(lambda v: f"{v:.4g}")
            tbl = ax.table(cellText=shown.values, colLabels=list(shown.columns), loc="upper left")
            tbl.auto_set_font_size(False)
            tbl.set_fontsize(7)
            pdf.savefig(fig)
            plt.close(fig)
            if figure is not None:
                pdf.savefig(figure)
                plt.close(figure)
    return path
