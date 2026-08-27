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

# shared collinear ray: pristine, m1, p1 all fall on this line (0.17 deg deviation)
pr = np.array(pristine); p1v = np.array(p1)
direction = p1v - pr
line_pts = np.array([pr + t*direction for t in np.linspace(-0.35, 1.15, 2)])
ax.plot(line_pts[:, 0], line_pts[:, 1], "--", color="0.5", lw=1, zorder=1,
        label="shared collinear ray (0.17° deviation)")

# confirmed single continuous valley: pristine -> m1 (chirality -1, compression)
xs = [pristine[0], m1[0]]
ys = [pristine[1], m1[1]]
ax.plot(xs, ys, "-", color="tab:blue", lw=2.5, zorder=2)
ax.annotate("", xy=m1, xytext=pristine,
            arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=2, shrinkA=8, shrinkB=8))

# basin boundary along the ray where pristine/m1's own attractor hands off
# to the (unrelated) k3=0.5 nodal plane, instead of continuing toward p1
boundary = pr + 0.7*direction
ax.scatter(*boundary, s=60, color="0.4", marker="x", zorder=4)
ax.annotate("basin boundary (t~0.7-0.8):\npristine/m1's attractor hands off\nto the k3=0.5 nodal plane here —\nNOT toward p1's point",
            xy=boundary, xytext=(0.03, 0.30), fontsize=8, color="0.35",
            arrowprops=dict(arrowstyle="->", color="0.4", lw=1))

ax.scatter(*pristine, s=140, color="tab:blue", zorder=3, marker="o",
           label=r"pristine (0%): $\chi=-1$")
ax.scatter(*m1, s=140, color="tab:blue", zorder=3, marker="s",
           label=r"m1 (nominal $-1\%$ biaxial): $\chi=-1$")
ax.scatter(*p1, s=140, color="tab:red", zorder=3, marker="^",
           label=r"p1 (nominal $+1\%$ biaxial): $\chi=+1$")

ax.annotate("confirmed single continuous valley\n(same point, same chirality)",
            xy=((pristine[0]+m1[0])/2, (pristine[1]+m1[1])/2),
            xytext=(0.17, 0.06), fontsize=9, color="tab:blue",
            arrowprops=dict(arrowstyle="->", color="tab:blue", lw=1))

ax.annotate("collinear with pristine/m1 (0.17°) but\nchirality flips (+1) and pristine/m1 have\nNO attractor here (basin goes to the plane\ninstead) — likely a separate touching event\nsharing the same 'soft' crystallographic\ndirection, not a continuous trajectory",
            xy=p1, xytext=(0.30, 0.15), fontsize=8.5, color="tab:red",
            arrowprops=dict(arrowstyle="->", color="tab:red", lw=1))

ax.annotate("m2.5 (verified $-2.5\\%$ biaxial,\na/c match to 0.0004%p):\nno Weyl point survives\n(band17-18 gap fully on\nk3=0.5 nodal plane)",
            xy=(0.5, 0.5), xytext=(0.34, 0.40), fontsize=8.5, color="0.35",
            ha="left")

ax.set_xlim(0, 0.5)
ax.set_ylim(0, 0.5)
ax.set_xlabel(r"$k_1$ (reduced)")
ax.set_ylabel(r"$k_3$ (reduced)")
ax.set_title("Band 17–18 gap Weyl points on the $k_2=0$ line\n(pristine / m1 / p1_mirror, non-transposed Born, LOTO_method='phonopy')")
ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig("results/weyl_trend/weyl_trend_k2eq0.png", dpi=170)
print("done")
