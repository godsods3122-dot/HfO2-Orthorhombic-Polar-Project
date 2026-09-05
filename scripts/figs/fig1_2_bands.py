#!/usr/bin/env python3
"""Fig 1: parent_pristine 전 경로 포논 분산.  Fig 2: band 17-18 확대 + 하이라이트."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY, BLUE
setup()
import matplotlib.pyplot as plt

NB, NK = 36, 120
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
d = np.loadtxt('figs/bulkek_parent_pristine.dat')
E = d[:, 1].reshape(NB, -1)
x = d[:len(d) // NB, 0]
nseg = len(LAB) - 1
ticks = [x[i * NK] for i in range(nseg)] + [x[-1]]

def frame(ax, ylim=None):
    for t in ticks[1:-1]:
        ax.axvline(t, color=GREY, lw=0.8, alpha=0.55)
    ax.set_xticks(ticks); ax.set_xticklabels(LAB)
    ax.set_xlim(x[0], x[-1])
    if ylim: ax.set_ylim(*ylim)
    ax.set_ylabel('Frequency (THz)')

# ---------------- Fig 1
fig, ax = plt.subplots(figsize=(9, 5.6))
for b in range(NB):
    ax.plot(x, E[b], color=GREY, lw=1.1)
frame(ax, (min(0, E.min() * 1.05), E.max() * 1.04))
ax.set_title('HfO$_2$ $Pca2_1$ — phonon dispersion (unstrained)')
fig.savefig('figs/fig1_dispersion_parent_pristine.png')
plt.close(fig)

# ---------------- Fig 2
g = E[17] - E[16]
i = int(np.argmin(g))
fig, ax = plt.subplots(figsize=(9, 5.6))
for b in range(NB):
    if b in (16, 17): continue
    ax.plot(x, E[b], color='#b9b9b9', lw=1.0)
ax.plot(x, E[16], color=BLUE, lw=2.4, label='band 17')
ax.plot(x, E[17], color='#63a4ff', lw=2.4, label='band 18')
lo = min(E[16].min(), E[17].min()) - 0.9
hi = max(E[16].max(), E[17].max()) + 0.9
frame(ax, (lo, hi))
ax.annotate('min gap on path\n%.3f THz' % g[i], xy=(x[i], (E[16][i] + E[17][i]) / 2),
            xytext=(x[i] + 0.32 * (x[-1] - x[0]) * 0.12, hi - 0.55),
            arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.6), color='#c0392b', fontsize=12)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2,
          frameon=False, fontsize=12)
ax.set_title('Bands 17 / 18 along the standard path')
fig.savefig('figs/fig2_bands17_18_highlight.png')
plt.close(fig)
print('fig1, fig2 저장.  경로 위 최소 gap = %.4e THz at %s' % (g[i], LAB))
print('band17 범위 %.3f~%.3f, band18 %.3f~%.3f THz' % (E[16].min(), E[16].max(), E[17].min(), E[17].max()))
