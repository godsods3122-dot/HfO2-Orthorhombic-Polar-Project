#!/usr/bin/env python3
"""발표용 단일 패널 — 노드를 지나는 밴드 17/18 절단 하나만, 폰트 크게.

fig2 의 (b)/(c) 를 각각 독립 그림으로 뽑는다. 축은 Δk 가 아니라 **절대 k 좌표**로
쓴다 (노드가 k_a=0.14649, k_b=0.07085 에 있으므로 눈금이 그 값을 지난다).
교점 레이블은 χ 만 남긴다.

사용: python3 fig2_cut_ppt.py
출력: figs/fig2_cut_ka.png   (k_a 방향 절단, k_b·k_c 고정)
      figs/fig2_cut_kb.png   (k_b 방향 절단, k_a·k_c 고정)
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import GREY, BLUE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NB = 36
NODE = (0.1464927, 0.0708493)      # (k_a, k_b)
WEYL_E = 10.086936
CHI = '+1'
HALF = 0.035                        # bulkek_nodecut.dat 이 노드 ±0.035 를 훑는다
RED, LIGHT = '#c0392b', '#67b0ff'

plt.rcParams.update({
    'font.size': 28, 'axes.labelsize': 36, 'axes.titlesize': 36,
    'xtick.labelsize': 30, 'ytick.labelsize': 30, 'legend.fontsize': 28,
    'axes.linewidth': 2.2, 'xtick.major.width': 2.2, 'ytick.major.width': 2.2,
    'xtick.major.size': 8, 'ytick.major.size': 8,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

c = np.loadtxt('figs/bulkek_nodecut.dat')
Ec = c[:, 1].reshape(NB, -1)
n = Ec.shape[1] // 2

for tag, sl, axis, fixed in (
        ('ka', slice(0, n),    'a', '$k_b$ = %.5f,  $k_c$ = 0' % NODE[1]),
        ('kb', slice(n, None), 'b', '$k_a$ = %.5f,  $k_c$ = 0' % NODE[0])):
    EE = Ec[:, sl]
    k0 = NODE[0] if axis == 'a' else NODE[1]
    k = np.linspace(k0 - HALF, k0 + HALF, EE.shape[1])
    i = int(np.argmin(EE[17] - EE[16]))

    fig, ax = plt.subplots(figsize=(12.4, 7.8))
    for b in range(NB):
        if b not in (16, 17):
            ax.plot(k, EE[b], color='#c8c8c8', lw=1.6)
    ax.axhline(WEYL_E, color=RED, lw=9.0, alpha=0.30, zorder=1)
    ax.plot(k, EE[16], color=BLUE, lw=4.2, label='band 17', zorder=3)
    ax.plot(k, EE[17], color=LIGHT, lw=4.2, label='band 18', zorder=3)
    ax.plot(k[i], EE[16][i], 'o', ms=20, mfc='none', mec=RED, mew=4.0, zorder=4)
    ax.annotate('$\\chi = %s$' % CHI, (k[i], EE[16][i]),
                textcoords='offset points', xytext=(0, -82), ha='center',
                fontsize=32, color=RED, fontweight='bold', zorder=5)

    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(EE[16].min() - 0.012, EE[17].max() + 0.012)
    ax.set_xlabel('$k_%s$  (reduced)' % axis, labelpad=10)
    ax.set_ylabel('Frequency (THz)', labelpad=10)
    # 고정 파라미터는 그림 위 제목 자리로 — 크게 키워도 밴드와 안 겹친다
    ax.set_title(fixed, fontsize=58, color='#555', pad=16)
    ax.legend(loc='upper left', frameon=False, handlelength=1.6)
    out = 'figs/fig2_cut_%s.png' % tag
    fig.savefig(out)
    plt.close(fig)
    print('wrote %s   node at k_%s = %.5f,  %.6f THz' % (out, axis, k[i], EE[16][i]))
