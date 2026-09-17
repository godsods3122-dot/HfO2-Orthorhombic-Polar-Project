#!/usr/bin/env python3
"""발표용 — 바일 4점 궤도의 strain 추이. 라벨은 strain 과 chirality 뿐.

미러 배치 세 점만 쓴다 (배치를 섞지 않는다). 좌표는 표준 (a, b, c), k_c = 0.
chirality 는 레포/Simphony 규약. 거울이 χ 를 뒤집으므로 1·3 사분면이 한 부호,
2·4 사분면이 반대 부호다.

출력: figs/fig8_trend_ppt.png
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# (라벨, (k_a, k_b), 1사분면 χ, 마커)
PTS = [('$-0.84$ %',  (0.1294182, 0.0637717), -1, 'o'),
       ('unstrained', (0.1610282, 0.0975234), -1, 'D'),
       ('$+1.16$ %',  (0.3182116, 0.2643381), +1, 's')]
POS, NEG, GREY = '#c0392b', '#1f5fd0', '#888888'

plt.rcParams.update({
    'font.size': 26, 'axes.labelsize': 40, 'legend.fontsize': 27,
    'xtick.labelsize': 30, 'ytick.labelsize': 30,
    'axes.linewidth': 2.4, 'xtick.major.width': 2.4, 'ytick.major.width': 2.4,
    'xtick.major.size': 9, 'ytick.major.size': 9,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

fig, ax = plt.subplots(figsize=(12.6, 10.6))
ax.axhline(0, color=GREY, lw=2.0, ls='--', zorder=1)
ax.axvline(0, color=GREY, lw=2.0, ls='--', zorder=1)

for sx in (1, -1):
    for sy in (1, -1):
        p = [(sx * a, sy * b) for _, (a, b), _, _ in PTS]
        for q, r in zip(p[:-1], p[1:]):
            ax.annotate('', xy=r, xytext=q,
                        arrowprops=dict(arrowstyle='-|>', color='#3a3a3a', lw=2.6,
                                        shrinkA=15, shrinkB=15, alpha=0.8), zorder=2)
        for (lab, (a, b), chi, m) in PTS:
            c = POS if chi * sx * sy > 0 else NEG
            ax.plot(sx * a, sy * b, m, ms=24, color=c, mec='k', mew=1.6, zorder=5)

ax.set_xlim(-0.44, 0.44); ax.set_ylim(-0.40, 0.40)
ax.set_aspect('equal')
ax.set_xlabel('$k_a$', labelpad=10); ax.set_ylabel('$k_b$', labelpad=6)
ax.set_xticks([-0.4, -0.2, 0.0, 0.2, 0.4])
ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4])

hs = [Line2D([], [], ls='', marker=m, ms=20, color='#c9c9c9', mec='k', mew=1.4, label=lab)
      for lab, _, _, m in PTS]
hc = [Line2D([], [], ls='', marker='o', ms=20, color=POS, mec='k', mew=1.4, label='$\\chi = +1$'),
      Line2D([], [], ls='', marker='o', ms=20, color=NEG, mec='k', mew=1.4, label='$\\chi = -1$')]
l1 = ax.legend(handles=hs, loc='upper center', bbox_to_anchor=(0.5, -0.115), ncol=3,
               frameon=False, handletextpad=0.35, columnspacing=1.3)
ax.add_artist(l1)
ax.legend(handles=hc, loc='upper center', bbox_to_anchor=(0.5, -0.215), ncol=2,
          frameon=False, handletextpad=0.5, columnspacing=3.0)

fig.savefig('figs/fig8_trend_ppt.png')
print('wrote figs/fig8_trend_ppt.png')
