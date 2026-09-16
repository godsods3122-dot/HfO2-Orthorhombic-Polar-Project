#!/usr/bin/env python3
"""발표용 — fig2 의 (a) 표준 경로 밴드 그림만, 폰트 크게 + 고대칭점 볼드.

fig2_node.py 의 (a) 와 같은 내용이다. 두 가지 주의는 그대로 적용된다:
  - 안쪽 Γ 는 LO-TO 방향 의존으로 불연속 (band 20 이 1.32 THz 튄다) → NaN 으로 끊는다
  - Γ-X 위 교점은 격자가 비껴가므로 phonopy 정밀화 값(t=0.185742)을 쓴다

출력: figs/fig2_path_ppt.png
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import GREY, BLUE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

NB, NK = 36, 120
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
GAMMA_BREAK = 4 * NK
WEYL_E = 10.086936
XC_FRAC, XC_E = 0.185742 / 0.5, 10.192507
RED, LIGHT = '#c0392b', '#67b0ff'

plt.rcParams.update({
    'font.size': 28, 'axes.labelsize': 38, 'legend.fontsize': 28,
    'ytick.labelsize': 32, 'xtick.labelsize': 36,
    'axes.linewidth': 2.4, 'xtick.major.width': 2.4, 'ytick.major.width': 2.4,
    'xtick.major.size': 9, 'ytick.major.size': 9,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'ytick.right': True,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

d = np.loadtxt('figs/bulkek_parent_pristine.dat')
E = d[:, 1].reshape(NB, -1); x = d[:len(d) // NB, 0]

# 안쪽 Γ 에서 **한 번만** 끊는다. NaN 만 넣으면 단차처럼 보이므로, 그 뒤 x 를 통째로
# 밀어 실제로 벌어진 흰 틈을 만든다 (데이터 값은 하나도 건드리지 않는다).
GAP = 0.022 * (x[-1] - x[0])
xs = x.copy(); xs[GAMMA_BREAK:] += GAP
ticks = [xs[i * NK] for i in range(len(LAB) - 1)] + [xs[-1]]
ticks[4] = x[GAMMA_BREAK] + GAP / 2            # Γ 눈금은 틈의 한가운데
xp = np.insert(xs, GAMMA_BREAK, np.nan)
Ep = np.insert(E, GAMMA_BREAK, np.nan, axis=1)

fig, ax = plt.subplots(figsize=(15.0, 9.2))
for b in range(NB):
    if b not in (16, 17):
        ax.plot(xp, Ep[b], color='#c8c8c8', lw=1.8)
ax.plot([ticks[0], ticks[1]], [WEYL_E, WEYL_E], color=RED, lw=11.0, alpha=0.30,
        solid_capstyle='butt', zorder=1)
ax.plot(xp, Ep[16], color=BLUE, lw=4.4, label='band 17', zorder=3)
ax.plot(xp, Ep[17], color=LIGHT, lw=4.4, label='band 18', zorder=3)
for n, t in enumerate(ticks[1:-1], start=1):
    if n == 4:                                  # 끊은 자리 — 틈 양쪽 경계에 긋는다
        ax.axvline(x[GAMMA_BREAK], color=GREY, lw=1.4, alpha=0.55)
        ax.axvline(x[GAMMA_BREAK] + GAP, color=GREY, lw=1.4, alpha=0.55)
    else:
        ax.axvline(t, color=GREY, lw=1.4, alpha=0.55)
ax.plot(ticks[0] + XC_FRAC * (ticks[1] - ticks[0]), XC_E, 'o',
        ms=20, mfc='none', mec=RED, mew=4.0, zorder=4)

ax.set_xticks(ticks)
ax.set_xticklabels(LAB, fontsize=42, fontweight='bold')   # 고대칭 경로 볼드
ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(9.1, 11.95)
ax.set_ylabel('Frequency (THz)', labelpad=12)
ax.tick_params(axis='x', pad=10)

h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], color=RED, lw=11.0, alpha=0.30),
      Line2D([], [], ls='none', marker='o', ms=16, mfc='none', mec=RED, mew=3.6)]
l += ['Weyl node  %.3f THz' % WEYL_E,
      'bands 17/18 meet  %.3f THz' % XC_E]
ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.135), ncol=2,
          frameon=False, handlelength=2.0, columnspacing=3.0, labelspacing=0.6)

fig.savefig('figs/fig2_path_ppt.png')
print('wrote figs/fig2_path_ppt.png')
