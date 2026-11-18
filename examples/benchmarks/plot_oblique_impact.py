"""Figure of the oblique impact benchmark: python plot_oblique_impact.py oblique_impact.csv

Copyright (c) 2026 Karol Tucki. MIT License.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLOURS = ["#0072B2", "#D55E00", "#009E73"]
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.45, "lines.linewidth": 0.6, "lines.markersize": 3.5,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.45, "ytick.major.width": 0.45,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.linewidth": 0.3, "grid.color": "#d0d0d0",
    "legend.frameon": False, "savefig.dpi": 600, "svg.fonttype": "none", "pdf.fonttype": 42,
})
src = sys.argv[1] if len(sys.argv) > 1 else "oblique_impact.csv"
d = np.genfromtxt(src, delimiter=",", names=True)
fig, (a, b) = plt.subplots(1, 2, figsize=(17.4 / 2.54, 2.6), gridspec_kw=dict(wspace=0.32))
a.plot(d["vt_vn"], d["wR_vn_theory"], color="black", lw=0.6, label="rigid-body theory")
k = slice(None, None, 2)
a.plot(d["vt_vn"][k], d["wR_vn_sim"][k], "o", mfc="white", mec=COLOURS[0], mew=0.6, label="TillBridge")
a.set_ylabel("rebound $\\omega R / v_n$ (–)"); a.set_ylim(0, 1.4); a.set_yticks(np.arange(0, 1.41, 0.2))
b.plot(d["vt_vn"], d["vt_vn"] - d["vt_out_vn_theory"], color="black", lw=0.6, label="rigid-body theory")
b.plot(d["vt_vn"][k], (d["vt_vn"] - d["vt_out_vn_sim"])[k], "s", mfc="white", mec=COLOURS[1], mew=0.6, label="TillBridge")
b.set_ylabel("tangential velocity loss $(v_t - v_t')/v_n$ (–)"); b.set_ylim(0, 0.6); b.set_yticks(np.arange(0, 0.61, 0.1))
for ax, l in ((a, "a"), (b, "b")):
    ax.set_xlim(0, 4); ax.set_xticks(np.arange(0, 4.01, 1)); ax.set_xlabel("incident $v_t / v_n$ (–)")
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    ax.text(0.5, -0.27, f"({l})", transform=ax.transAxes, ha="center", va="top")
out = src.rsplit(".", 1)[0]
for ext in ("png", "svg"):
    fig.savefig(f"{out}.{ext}", bbox_inches="tight")
print("written", out + ".png")
