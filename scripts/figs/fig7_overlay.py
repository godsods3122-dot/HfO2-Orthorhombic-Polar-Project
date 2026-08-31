#!/usr/bin/env python3
"""Fig 7: biaxial strain 에 따른 밴드 추이 겹침 (rainbow 그라데이션, 반투명)."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

NB, NK = 36, 120
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
SETS = [('parent_pristine', 0.0, 'pristine  (0 %)'),
        ('m08', -0.8, '−0.8 %'),
        ('m10', -1.0, '−1.0 %'),
        ('m25', -2.5, '−2.5 %'),
        ('m30', -3.0, '−3.0 %')]
cmap = plt.get_cmap('rainbow')
norm = Normalize(vmin=-3.0, vmax=0.0)

fig, axs = plt.subplots(1, 2, figsize=(15.4, 6.0), gridspec_kw=dict(width_ratios=[1.55, 1], wspace=0.2))
x = None
for tag, eps, lab in SETS:
    d = np.loadtxt('figs/bulkek_%s.dat' % tag)
    E = d[:, 1].reshape(NB, -1)
    if x is None:
        x = d[:len(d) // NB, 0]
        ticks = [x[i * NK] for i in range(len(LAB) - 1)] + [x[-1]]
    c = cmap(norm(eps))
    for b in range(NB):
        axs[0].plot(x, E[b], color=c, lw=1.15, alpha=0.62, zorder=2)
    axs[0].plot([], [], color=c, lw=2.6, label=lab)
    for b in (16, 17):
        axs[1].plot(x, E[b], color=c, lw=2.0, alpha=0.85)
    axs[1].plot([], [], color=c, lw=2.6, label=lab)

for ax, ylim, ttl in ((axs[0], (0, 24.4), '(a)  all 36 branches'),
                      (axs[1], (9.0, 12.4), '(b)  bands 17 / 18')):
    for t in ticks[1:-1]:
        ax.axvline(t, color=GREY, lw=0.8, alpha=0.45, zorder=1)
    ax.set_xticks(ticks); ax.set_xticklabels(LAB); ax.set_xlim(x[0], x[-1])
    ax.set_ylim(*ylim); ax.set_ylabel('Frequency (THz)')
    ax.set_title(ttl, loc='left')
axs[0].legend(loc='upper left', ncol=2, framealpha=0.94, fontsize=11.5)
axs[1].legend(loc='lower left', framealpha=0.94, fontsize=11.5)
sm = ScalarMappable(norm=norm, cmap=cmap); sm.set_array([])
cb = fig.colorbar(sm, ax=axs, fraction=0.022, pad=0.015)
cb.set_label('biaxial strain (%)')
fig.suptitle('HfO$_2$ $Pca2_1$ — phonon bands under biaxial compression', y=0.98, fontsize=15)
fig.savefig('figs/fig7_strain_overlay.png')
print('fig7 저장')
