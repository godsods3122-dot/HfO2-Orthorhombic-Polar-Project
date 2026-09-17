#!/usr/bin/env python3
"""발표용 단일 패널 — 미러 −1 % 노드를 지나는 밴드 17/18 절단.

데이터는 gen_m1_data.py 의 figs/m1_figdata.npz (phonopy 직접 계산).
축은 Δk 가 아니라 절대 k 좌표라 눈금이 노드 값을 지나간다.
교점 레이블은 χ 만 — 부호는 레포/Simphony 규약이다 (자체 Berry flux 는 전체
부호가 반대이므로 상대 부호만 비교할 것; chi_robust.py 참조).

출력: figs/fig2_cut_ka.png, figs/fig2_cut_kb.png

⚠️ chirality 부호 — **parent 기준으로 표기**한다.
계산은 미러 구조로 했고 (Simphony WeylChirality_calc, runs/<src>/band17/PN.out),
미러와 parent 는 det = −1 (improper) 변환으로 연결되므로 chirality 가 뒤집힌다.
따라서 미러 계산값의 부호를 반전해 적는다.

  구조              native (k1,k2,k3)          Simphony(미러)   표기(parent)
  pristine_mirror   (0.09752, 0, 0.16103)          −1              +1
  m1_mirror         (0.06377, 0, 0.12942)          −1              +1
  p1_mirror         (0.26434, 0, 0.31821)          +1              −1

색 규약: 양수 빨강 / 음수 파랑.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import BLUE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BAND = 17
CHI = '+1'                      # (+,+) 노드, parent 기준 (미러 −1 의 반전)
RED, LIGHT = '#c0392b', '#67b0ff'

d = np.load('figs/m1_figdata.npz', allow_pickle=True)
NODE = d['node']
WEYL_E = float(d['cone_E17'].ravel()[np.argmin(d['cone_E18'] - d['cone_E17'])])

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

for name, other in (('a', ('b', 1)), ('b', ('a', 0))):
    k = d['cut_k_' + name]
    EE = d['cut_E_' + name].T
    fixed = '$k_%s$ = %.5f,  $k_c$ = 0' % (other[0], NODE[other[1]])
    i = int(np.argmin(EE[BAND] - EE[BAND-1]))

    fig, ax = plt.subplots(figsize=(12.4, 7.8))
    for b in range(EE.shape[0]):
        if b not in (BAND-1, BAND):
            ax.plot(k, EE[b], color='#c8c8c8', lw=1.6)
    ax.axhline(WEYL_E, color=RED, lw=9.0, alpha=0.30, zorder=1)
    ax.plot(k, EE[BAND-1], color=BLUE, lw=4.2, label='band 17', zorder=3)
    ax.plot(k, EE[BAND], color=LIGHT, lw=4.2, label='band 18', zorder=3)
    ax.plot(k[i], EE[BAND-1][i], 'o', ms=20, mfc='none', mec=RED, mew=4.0, zorder=4)
    ax.annotate('$\\chi = %s$' % CHI, (k[i], EE[BAND-1][i]),
                textcoords='offset points', xytext=(0, -82), ha='center',
                fontsize=32, color=RED, fontweight='bold', zorder=5)

    ax.set_xlim(k[0], k[-1])
    ax.set_ylim(EE[BAND-1].min() - 0.012, EE[BAND].max() + 0.012)
    ax.set_xlabel('$k_%s$  (reduced)' % name, labelpad=10)
    ax.set_ylabel('Frequency (THz)', labelpad=10)
    ax.set_title(fixed, fontsize=58, color='#555', pad=16)
    ax.legend(loc='upper left', frameon=False, handlelength=1.6)
    out = 'figs/fig2_cut_k%s.png' % name
    fig.savefig(out); plt.close(fig)
    print('wrote %s   node k_%s = %.5f,  %.6f THz,  gap %.2e'
          % (out, name, k[i], EE[BAND-1][i], (EE[BAND] - EE[BAND-1])[i]))
