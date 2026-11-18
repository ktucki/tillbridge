"""Figures of the chisel example: python make_figures.py [OUTDIR]  (run chisel_example first).

Reads only the text and CSV files written by TillBridge; needs NumPy and Matplotlib.
Writes PNG (600 dpi) and SVG files to OUTDIR.

Copyright (c) 2026 Karol Tucki. MIT License.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.ticker import MaxNLocator

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "out")
COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]   # Okabe-Ito
X0, X1 = 0.30, 0.75
W1, W2 = 8.4 / 2.54, 17.4 / 2.54   # single and double column width (in)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.linewidth": 0.45, "lines.linewidth": 0.6, "lines.markersize": 3.5,
    "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.45, "ytick.major.width": 0.45,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.linewidth": 0.3, "grid.color": "#d0d0d0",
    "legend.frameon": False, "savefig.dpi": 600, "axes.prop_cycle": matplotlib.cycler(color=COLOURS),
    "svg.fonttype": "none", "pdf.fonttype": 42,
})


def nice(ax, axis="y", zero=False, nbins=6):
    """Axis from a labelled tick to a labelled tick."""
    lo, hi = ax.get_ylim() if axis == "y" else ax.get_xlim()
    if zero:
        lo = min(0.0, lo)
    t = MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10]).tick_values(lo, hi)
    step = t[1] - t[0]
    a, b = np.floor(lo / step + 1e-9) * step, np.ceil(hi / step - 1e-9) * step
    ticks = np.arange(a, b + step / 2, step)
    (ax.set_ylim if axis == "y" else ax.set_xlim)(a, b)
    (ax.set_yticks if axis == "y" else ax.set_xticks)(ticks)


def save(fig, name):
    for ext in ("png", "svg"):
        fig.savefig(OUT / f"{name}.{ext}", bbox_inches="tight")
    plt.close(fig)


def keyvalue(path):
    d = {}
    for line in open(path):
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                d[k.strip()] = float(v)
            except ValueError:
                d[k.strip()] = v.strip()
    return d


def force(prefix):
    return np.loadtxt(OUT / f"{prefix}_force.csv", delimiter=",", skiprows=1)


def tool(path):
    return np.loadtxt(path, skiprows=1).reshape(-1, 3, 3)


def calibration_log(path):
    rows = []
    for line in open(path):
        if "|" not in line:
            continue
        a, b = line.split("|")
        f = a.split()
        r = dict(kv.split("=") for kv in b.split())
        rows.append(dict(label=f[0], src=f[1], x=[float(v) for v in f[2:]], **{k: float(v) for k, v in r.items()}))
    return rows


S = keyvalue(OUT / "summary.txt")
CAL = calibration_log(OUT / "calibration.log")

# ------------------------------------------------------------ Fig. 2 calibration
x = np.array([c["x"] for c in CAL]); lhs = np.array([c["src"] == "lhs" for c in CAL])
D = np.array([c["draft"] for c in CAL]); V = np.array([c["vertical"] for c in CAL])
fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 2.7), gridspec_kw=dict(wspace=0.32))
a.plot(x[lhs, 0] / 1e3, x[lhs, 1], "o", mfc="white", mec=COLOURS[0], mew=0.6, label="start design")
a.plot(x[~lhs, 0] / 1e3, x[~lhs, 1], "s", mfc="white", mec=COLOURS[1], mew=0.6, label="Bayesian steps")
a.plot(S["k_coh_true"] / 1e3, S["mu_pw_true"], "x", color="black", ms=6, mew=0.8, label="reference")
a.plot(S["k_coh_estimate"] / 1e3, S["mu_pw_estimate"], "+", color=COLOURS[2], ms=9, mew=0.8, label="estimate")
a.errorbar(S["k_coh_calibrated"] / 1e3, S["mu_pw_measured"], xerr=[[(S["k_coh_calibrated"] - S["k_coh_low"]) / 1e3], [(S["k_coh_high"] - S["k_coh_calibrated"]) / 1e3]],
           fmt="D", ms=3.5, mfc="white", mec=COLOURS[3], ecolor=COLOURS[3], elinewidth=0.6, capsize=2, mew=0.6, label="$\\mu_{pw}$ measured")
mp = np.genfromtxt(OUT / "calibration_map.csv", delimiter=",", names=True)
nc, nm = len(np.unique(mp["k_coh"])), len(np.unique(mp["mu_pw"]))
a.contour(mp["k_coh"].reshape(nc, nm) / 1e3, mp["mu_pw"].reshape(nc, nm), mp["chi2"].reshape(nc, nm), levels=[2.30], colors="black", linewidths=0.6)
a.contour(mp["k_coh"].reshape(nc, nm) / 1e3, mp["mu_pw"].reshape(nc, nm), mp["chi2"].reshape(nc, nm), levels=[6.18], colors="black", linewidths=0.6, linestyles=":")
a.set_xlim(0, 180); a.set_xticks(np.arange(0, 181, 30)); a.set_ylim(0.2, 0.9); a.set_yticks(np.arange(0.2, 0.91, 0.1))
a.set_xlabel("cohesion $k_{coh}$ (kPa)"); a.set_ylabel("soil–steel friction $\\mu_{pw}$ (–)")
a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.3, columnspacing=1.0)
a.grid(False)
n = np.arange(1, len(CAL) + 1)
b.plot(n, D, "o-", mfc="white", mew=0.6, label="draft force")
b.plot(n, V, "s-", mfc="white", mew=0.6, label="vertical force")
b.axhline(S["draft_reference"], color=COLOURS[0], ls=":", lw=0.6, label="reference values")
b.axhline(S["vertical_reference"], color=COLOURS[1], ls=":", lw=0.6)
b.set_xlabel("run (–)"); b.set_ylabel("force (N)")
b.set_xlim(0, 2 * int(np.ceil((len(CAL) + 1) / 2))); b.set_xticks(np.arange(0, b.get_xlim()[1] + 1, 2))
nice(b, "y", zero=True)
b.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.3, columnspacing=1.0)
for ax, l in ((a, "a"), (b, "b")):
    ax.text(0.5, -0.27, f"({l})", transform=ax.transAxes, ha="center", va="top")
save(fig, "Fig2_calibration")

# ------------------------------------------------------------ Fig. 3 force history and wear index
fb = force("calibrated")
T = tool(OUT / "tine.txt")
area = 0.5 * np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1)
Pt = np.loadtxt(OUT / "calibrated_triW.csv", delimiter=",", skiprows=1)[:, 1]
qf = Pt / area / 1e3   # kW/m2
fig, (a, b) = plt.subplots(1, 2, figsize=(W2, 2.7), gridspec_kw=dict(wspace=0.32, width_ratios=[1.35, 1]))
m = fb[:, 0] >= 0.1
a.plot(fb[m, 5], -fb[m, 1], color=COLOURS[0], lw=0.45, label="draft")
a.plot(fb[m, 5], -fb[m, 3], color=COLOURS[1], lw=0.45, label="vertical force")
a.plot([X0, X1], [S["draft_reference"]] * 2, color="black", ls=":", lw=0.6, label="reference means")
a.plot([X0, X1], [S["vertical_reference"]] * 2, color="black", ls=":", lw=0.6)
a.set_xlim(0, 0.8); a.set_xticks(np.arange(0, 0.81, 0.2))
a.set_xlabel("tool displacement (m)"); a.set_ylabel("force (N)")
nice(a, "y")
a.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, handletextpad=0.3, columnspacing=1.0)
# front face of the tine seen along the travel direction (triangles facing +x)
nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); front = nrm[:, 0] / np.linalg.norm(nrm, axis=1) > 0.3
yz = T[front][:, :, [1, 2]] * 100
pc = PolyCollection(yz, array=qf[front], cmap="viridis", edgecolors="face", linewidths=0.2)
b.add_collection(pc); b.set_aspect("equal"); b.grid(False)
b.set_xlim(18, 32); b.set_xticks(np.arange(18, 32.1, 2)); b.set_ylim(8, 30); b.set_yticks(np.arange(8, 30.1, 2))
b.set_xlabel("y (cm)"); b.set_ylabel("z (cm)")
cb = fig.colorbar(pc, ax=b, fraction=0.08, pad=0.04)
cb.set_label("frictional power density (kW/m$^2$)", fontsize=8); cb.ax.tick_params(labelsize=7, width=0.45); cb.outline.set_linewidth(0.45)
b.axhline(20, color="black", lw=0.45, ls="--")
b.text(b.get_xlim()[1], 20.4, "soil surface", ha="right", va="bottom", fontsize=7)
for ax, l in ((a, "a"), (b, "b")):
    ax.text(0.5, -0.27, f"({l})", transform=ax.transAxes, ha="center", va="top")
save(fig, "Fig3_force_wear")

# ------------------------------------------------------------ Fig. 4 shank stress
lines = open(OUT / "shank.vtk").read().split("\n")
i = next(k for k, l in enumerate(lines) if l.startswith("POINTS")); nn = int(lines[i].split()[1])
X = np.array([[float(v) for v in l.split()] for l in lines[i + 1:i + 1 + nn]])
j = next(k for k, l in enumerate(lines) if l.startswith("SCALARS von_mises"))
vm = np.array([float(l) for l in lines[j + 2:j + 2 + nn]])
k = next(k for k, l in enumerate(lines) if l.startswith("VECTORS displacement"))
U = np.array([[float(v) for v in l.split()] for l in lines[k + 1:k + 1 + nn]])
side = np.abs(X[:, 1] - X[:, 1].min()) < 1e-9   # side face y = -25 mm
fig, (a, b) = plt.subplots(1, 2, figsize=(W1 * 1.25, 3.0), gridspec_kw=dict(wspace=0.6))
tc = a.tricontourf(X[side, 0], X[side, 2], vm[side], levels=12, cmap="viridis")
a.set_aspect("equal"); a.grid(False); a.set_xlim(0, 20); a.set_xticks([0, 20]); a.set_ylim(0, 350); a.set_yticks(np.arange(0, 351, 50))
a.set_xlabel("x (mm)"); a.set_ylabel("distance from the tip (mm)")
cb = fig.colorbar(tc, ax=a, fraction=0.5, pad=0.15)
cb.set_label("von Mises stress (MPa)", fontsize=8); cb.ax.tick_params(labelsize=7, width=0.45); cb.outline.set_linewidth(0.45)
axis = (X[:, 0] - 10) ** 2 + X[:, 1] ** 2 < 1e-6
o = np.argsort(X[axis, 2])
b.plot(np.linalg.norm(U[axis], axis=1)[o], X[axis, 2][o], color=COLOURS[0])
b.set_ylim(0, 350); b.set_yticks(np.arange(0, 351, 50)); b.set_xlabel("|u| (mm)"); b.grid(True, axis="x"); b.grid(False, axis="y")
nice(b, "x", zero=True, nbins=3)
for ax, l in ((a, "a"), (b, "b")):
    ax.text(0.5, -0.13, f"({l})", transform=ax.transAxes, ha="center", va="top")
save(fig, "Fig4_shank")
print("figures written to", OUT)
