#!/usr/bin/env python3
"""발표용 — 미러 −1 % 의 Weyl 콘 3D, 폰트 크게, 레이블은 χ 만.

데이터는 gen_m1_data.py 의 figs/m1_figdata.npz (phonopy 직접 계산).
parent_pristine 판본(fig5_cone.py)과 달리 미러 계열이다.

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D

CHI = '+1'                      # (+,+) 노드, parent 기준 (미러 −1 의 반전)
BLU, ORG, RED = '#2b6cb0', '#dd6b20', '#c0392b'

d = np.load('figs/m1_figdata.npz', allow_pickle=True)
a, E17, E18 = d['cone_a'], d['cone_E17'], d['cone_E18']
N = len(a)
A, B = np.meshgrid(a, a, indexing='ij')
i = np.unravel_index(np.argmin(E18 - E17), E17.shape)
w0 = E17[i]; j = i[0]
print('미러 −1%% 콘  gap %.3e THz  E %.6f THz  node %s'
      % ((E18 - E17)[i], w0, d['node']))

plt.rcParams.update({
    'font.size': 22, 'axes.labelsize': 28, 'legend.fontsize': 21,
    'xtick.labelsize': 19, 'ytick.labelsize': 19,
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})
fig = plt.figure(figsize=(11.0, 9.4))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(A, B, E18, color=ORG, alpha=0.62, linewidth=0, antialiased=True,
                rstride=2, cstride=2, shade=True)
ax.plot_surface(A, B, E17, color=BLU, alpha=0.62, linewidth=0, antialiased=True,
                rstride=2, cstride=2, shade=True)
for E, c in ((E17, BLU), (E18, ORG)):
    ax.plot(a, np.full(N, 0.0), E[:, j], color=c, lw=4.6, zorder=9)
    ax.plot(np.full(N, 0.0), a, E[j, :], color=c, lw=4.6, ls='--', zorder=9)
ax.scatter([0], [0], [w0], s=190, c='k', depthshade=False, zorder=12)
ax.text2D(0.44, 0.75, '$\\chi = %s$' % CHI, transform=ax.transAxes,
          fontsize=40, color=RED, fontweight='bold', ha='right', va='center')
ax.annotate('', xy=(0.512, 0.512), xytext=(0.452, 0.730), xycoords='axes fraction',
            textcoords='axes fraction', arrowprops=dict(arrowstyle='-', color='#444', lw=2.0))
ax.set_xlabel('$\\Delta k_a$', labelpad=40)
ax.set_ylabel('$\\Delta k_b$', labelpad=40)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('Frequency (THz)', labelpad=46, rotation=90)
ax.view_init(elev=18, azim=-61)
ax.set_box_aspect((1, 1, 0.78), zoom=1.02)
ax.set_zlim(E17.min(), E18.max())
ax.tick_params(labelsize=19, pad=3); ax.tick_params(axis='z', pad=16)
ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(4))
ax.zaxis.set_major_locator(MaxNLocator(5))
ax.legend(handles=[Line2D([], [], color=BLU, lw=6, label='band 17'),
                   Line2D([], [], color=ORG, lw=6, label='band 18'),
                   Line2D([], [], color='k', lw=3.2, ls='-', label='cut along $\\Delta k_a$'),
                   Line2D([], [], color='k', lw=3.2, ls='--', label='cut along $\\Delta k_b$')],
          loc='upper center', bbox_to_anchor=(0.5, -0.045), ncol=2, frameon=False)
fig.subplots_adjust(left=0.0, right=0.92, bottom=0.14, top=1.02)
fig.savefig('figs/fig5a_weyl_cone.png', dpi=300, bbox_inches='tight')
print('wrote figs/fig5a_weyl_cone.png')
