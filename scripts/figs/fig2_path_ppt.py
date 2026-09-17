#!/usr/bin/env python3
"""발표용 — 미러 −1 % 구조의 표준 경로 밴드 (band 17/18 강조), 폰트 크게.

데이터는 gen_m1_data.py 가 phonopy 로 직접 만든 figs/m1_figdata.npz 다.
parent_pristine(일반 구조)과 섞지 않는다 — 두 계열은 주파수가 0.14 THz 어긋난다.

안쪽 Γ 는 LO-TO 비해석항의 방향 의존 때문에 값이 불연속이므로 NaN 하나로 끊는다.
Γ-X 는 거울면 위 선이라 band 17/18 이 만날 수 있고, 실제로 만난다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import GREY, BLUE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RED, LIGHT = '#c0392b', '#67b0ff'
d = np.load('figs/m1_figdata.npz', allow_pickle=True)
x, E = d['path_x'], d['path_E']
ticks = d['ticks']
LAB = [('Γ' if str(s) == 'G' else str(s)) for s in d['labels']]
GB = int(d['gamma_break'])
WEYL_E = float(d['cone_E17'].ravel()[np.argmin(d['cone_E18'] - d['cone_E17'])])
XC_T, XC_E = float(d['xc_t']), float(d['xc_E'])

plt.rcParams.update({
    'font.size': 28, 'axes.labelsize': 48, 'legend.fontsize': 28,
    'ytick.labelsize': 34, 'xtick.labelsize': 60,
    'axes.linewidth': 2.4, 'xtick.major.width': 2.4, 'ytick.major.width': 2.4,
    'xtick.major.size': 9, 'ytick.major.size': 9,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'ytick.right': True,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

xp = np.insert(x, GB, np.nan)
Ep = np.insert(E, GB, np.nan, axis=1)

fig, ax = plt.subplots(figsize=(15.0, 9.2))
for b in range(E.shape[0]):
    if b not in (16, 17):
        ax.plot(xp, Ep[b], color='#c8c8c8', lw=1.8)
ax.plot([ticks[0], ticks[1]], [WEYL_E, WEYL_E], color=RED, lw=11.0, alpha=0.30,
        solid_capstyle='butt', zorder=1)
ax.plot(xp, Ep[16], color=BLUE, lw=4.4, label='band 17', zorder=3)
ax.plot(xp, Ep[17], color=LIGHT, lw=4.4, label='band 18', zorder=3)
for t in ticks[1:-1]:
    ax.axvline(t, color=GREY, lw=1.4, alpha=0.55)
ax.plot(ticks[0] + (XC_T / 0.5) * (ticks[1] - ticks[0]), XC_E, 'o',
        ms=20, mfc='none', mec=RED, mew=4.0, zorder=4)

ax.set_xticks(ticks)
ax.set_xticklabels(LAB, fontsize=60)
ax.set_xlim(x[0], x[-1])
lo = min(Ep[16][~np.isnan(Ep[16])].min(), WEYL_E) - 0.45
hi = max(Ep[17][~np.isnan(Ep[17])].max(), XC_E) + 0.35
ax.set_ylim(lo, hi)
ax.set_ylabel('Frequency (THz)', labelpad=12)
ax.tick_params(axis='x', pad=10)

h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], color=RED, lw=11.0, alpha=0.30),
      Line2D([], [], ls='none', marker='o', ms=16, mfc='none', mec=RED, mew=3.6)]
l += ['Weyl node  %.3f THz' % WEYL_E, 'bands 17/18 meet  %.3f THz' % XC_E]
ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.135), ncol=2,
          frameon=False, handlelength=2.0, columnspacing=3.0, labelspacing=0.6)

fig.savefig('figs/fig2_path_ppt.png')
print('wrote figs/fig2_path_ppt.png   (미러 −1%%)  Weyl %.6f THz,  Γ-X 교점 %.6f THz'
      % (WEYL_E, XC_E))
