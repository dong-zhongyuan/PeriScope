"""Shared figure style (adapted from the lxh kidney project): sans-serif 6-7 pt text, 8 pt bold panel labels, 0.5 pt
axes, italic P, colour-blind-safe Okabe-Ito colours. Springer/ESR widths: 84 mm single, 174 mm double column.
Vector PDF + 300 dpi PNG."""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

MM = 1 / 25.4
W1, W2 = 84 * MM, 174 * MM
INK, MUTED, FAINT, GRID = "#222222", "#666666", "#BBBBBB", "#EBEBEB"
# Okabe-Ito
BLUE, ORANGE, GREEN, VERMIL, SKY, PURPLE, YELLOW = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7", "#F0E442"
BENIGN, MALIGNANT = BLUE, VERMIL
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6.5, "axes.titlesize": 7, "axes.titlelocation": "left", "axes.titlepad": 3.5,
    "axes.labelsize": 6.5, "axes.labelpad": 2, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "legend.frameon": False, "legend.handlelength": 1.4, "legend.borderaxespad": 0.3,
    "axes.linewidth": 0.5, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "lines.linewidth": 0.9, "lines.markersize": 3, "patch.linewidth": 0.5, "errorbar.capsize": 0,
    "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 300, "mathtext.default": "regular",
})


def label(ax, s, dx=-18, dy=4):
    """Panel label 'a' in 8 pt bold, dx, dy points from the axes' top-left corner."""
    return ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                       fontsize=8, fontweight="bold", ha="left", va="bottom", annotation_clip=False)


def save(fig, outdir, name):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{name}.pdf")
    fig.savefig(outdir / f"{name}.png")
    plt.close(fig)
    print("WROTE", name, flush=True)


def pfmt(p):
    P = "$\\it{P}$"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{P} = {m} × 10$^{{{int(e)}}}$"
    return f"{P} = {p:.2g}" if p < 0.1 else f"{P} = {p:.2f}"


def roc(y, s):
    """ROC points (fpr, tpr) with ties handled as one step."""
    y = np.asarray(y, bool)
    s = np.asarray(s, float)
    o = np.argsort(-s, kind="stable")
    y, s = y[o], s[o]
    tp, fp = np.cumsum(y), np.cumsum(~y)
    keep = np.r_[np.diff(s) != 0, True]
    return np.r_[0, fp[keep] / (~y).sum()], np.r_[0, tp[keep] / y.sum()]


def spread(vals, gap):
    """Label positions for values sorted in descending order, pushed apart to at least `gap` keeping their order."""
    y = np.array(vals, float)
    for _ in range(500):
        moved = False
        for i in range(len(y) - 1):
            d = y[i] - y[i + 1]
            if d < gap - 1e-12:
                s = (gap - d) / 2
                y[i] += s
                y[i + 1] -= s
                moved = True
        if not moved:
            break
    return list(y)
