#!/usr/bin/env python3
"""Fig 2: band 17/18 하이라이트 + 실제 Weyl 노드 통과 선."""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY, BLUE
setup()
import matplotlib.pyplot as plt

NB, NK = 36, 120
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
d = np.loadtxt('figs/bulkek_parent_pristine.dat')
E = d[:, 1].reshape(NB, -1); x = d[:len(d) // NB, 0]
ticks = [x[i * NK] for i in range(len(LAB) - 1)] + [x[-1]]

c = np.loadtxt('figs/bulkek_nodecut.dat')
Ec = c[:, 1].reshape(NB, -1); xc = c[:len(c) // NB, 0]
n = Ec.shape[1] // 2
xa, Ea = xc[:n], Ec[:, :n]
xb, Eb = xc[n:] - xc[n], Ec[:, n:]

fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2),
                         gridspec_kw=dict(width_ratios=[2.15, 1, 1], wspace=0.28))

ax = axes[0]
for b in range(NB):
    if b not in (16, 17):
        ax.plot(x, E[b], color='#c8c8c8', lw=1.0)
ax.plot(x, E[16], color=BLUE, lw=2.5, label='band 17')
ax.plot(x, E[17], color='#67b0ff', lw=2.5, label='band 18')
for t in ticks[1:-1]:
    ax.axvline(t, color=GREY, lw=0.8, alpha=0.5)
ax.axvspan(ticks[5], ticks[9], color='#ffd9d9', alpha=0.55, zorder=0)
ax.text((ticks[5] + ticks[9]) / 2, 11.62, '$k_c=\\frac{1}{2}$ nodal plane\n(2$_1$+T enforced, $\\chi$=0)',
        ha='center', va='top', fontsize=11, color='#a33')
ax.set_xticks(ticks); ax.set_xticklabels(LAB); ax.set_xlim(x[0], x[-1])
ax.set_ylim(9.1, 11.95); ax.set_ylabel('Frequency (THz)')
ax.legend(loc='lower left', frameon=True, framealpha=0.95)
ax.set_title('(a)  bands 17 / 18 on the standard path', loc='left')

for ax, xx, EE, ttl, sub in ((axes[1], xa, Ea, '(b)  cut along $k_a$', '$k_b$=0.0708, $k_c$=0'),
                             (axes[2], xb, Eb, '(c)  cut along $k_b$', '$k_a$=0.1465, $k_c$=0')):
    for b in range(NB):
        if b not in (16, 17):
            ax.plot(xx, EE[b], color='#c8c8c8', lw=1.0)
    ax.plot(xx, EE[16], color=BLUE, lw=2.6)
    ax.plot(xx, EE[17], color='#67b0ff', lw=2.6)
    g = EE[17] - EE[16]; i = int(np.argmin(g))
    ax.plot(xx[i], EE[16][i], 'o', ms=9, mfc='none', mec='#c0392b', mew=2.2)
    ax.set_xticks([xx[0], xx[i], xx[-1]]); ax.set_xticklabels(['−0.035', '0', '+0.035'])
    ax.set_xlabel('Δk (reduced)')
    ax.set_ylim(EE[16].min() - 0.09, EE[17].max() + 0.09)
    ax.set_title(ttl, loc='left'); ax.text(0.5, 0.03, sub, transform=ax.transAxes,
                                           ha='center', fontsize=11, color='#555')
axes[1].set_ylabel('Frequency (THz)')
fig.text(0.655, 0.955, 'Weyl node  $k$=(0.1465, 0.0708, 0),  10.087 THz = 41.72 meV,  $\\chi=+1$',
         ha='center', fontsize=12.5, color='#c0392b')
fig.savefig('figs/fig2_bands17_18_node.png')
print('fig2 저장. 노드 선 위 최소 gap = %.2e THz' % (Ec[17] - Ec[16]).min())
