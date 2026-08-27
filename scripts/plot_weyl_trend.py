#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(8, 8))

# band17-18 gap, k2=0 line, (k1,k3) — chirality -1 cluster (pristine -> m1, continuous)
pristine = (0.0975233, 0.1610283)
m1 = (0.0637718, 0.1294180)
# band17-18 gap, k2=0 line — chirality +1 point at p1 (verified NOT continuous with the above)
p1 = (0.2643385, 0.3182119)

# mirror-plane guides (chirality forced to 0)
for v in (0.0, 0.5):
    ax.axvline(v, color="0.75", lw=1, zorder=0)
    ax.axhline(v, color="0.75", lw=1, zorder=0)
ax.text(0.008, 0.015, "k1 mirror plane", color="0.55", fontsize=8, rotation=90, va="bottom")
ax.text(0.485, 0.008, "k3 mirror plane", color="0.55", fontsize=8, ha="right")

# continuous trend: pristine -> m1 (chirality -1, compression)
xs = [pristine[0], m1[0]]
ys = [pristine[1], m1[1]]
ax.plot(xs, ys, "-", color="tab:blue", lw=2, zorder=2)
ax.annotate("", xy=m1, xytext=pristine,
            arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=2, shrinkA=8, shrinkB=8))

ax.scatter(*pristine, s=140, color="tab:blue", zorder=3, marker="o",
           label=r"pristine (0%): $\chi=-1$")
ax.scatter(*m1, s=140, color="tab:blue", zorder=3, marker="s",
           label=r"m1 (nominal $-1\%$ biaxial): $\chi=-1$")
ax.scatter(*p1, s=140, color="tab:red", zorder=3, marker="^",
           label=r"p1 (nominal $+1\%$ biaxial): $\chi=+1$")

ax.annotate("continuous compression trend\n(same point, same chirality)",
            xy=((pristine[0]+m1[0])/2, (pristine[1]+m1[1])/2),
            xytext=(0.16, 0.09), fontsize=9, color="tab:blue",
            arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1))

ax.annotate("separate point — chirality flips (+1),\nverified NOT a continuation of the\npristine/m1 pair (direct gap evaluation)",
            xy=p1, xytext=(0.30, 0.22), fontsize=9, color="tab:red",
            arrowprops=dict(arrowstyle="->", color="tab:red", lw=1))

ax.annotate("m2.5 (verified $-2.5\\%$ biaxial,\na/c match to 0.0004%p):\nno Weyl point survives\n(band17-18 gap fully on\nk3=0.5 nodal plane)",
            xy=(0.5, 0.5), xytext=(0.34, 0.40), fontsize=8.5, color="0.35",
            ha="left")

ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)
ax.set_xlabel(r"$k_1$ (reduced)")
ax.set_ylabel(r"$k_3$ (reduced)")
ax.set_title("Band 17–18 gap Weyl points on the $k_2=0$ line\n(pristine / m1 / p1_mirror, non-transposed Born, LOTO_method='phonopy')")
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig("/tmp/claude-0/-home-user-HfO2-Orthorhombic-Polar-Project/c2573bf3-53de-56be-b0b9-1d949d54780e/scratchpad/weyl_trend_k2eq0.png", dpi=170)
print("done")
