#!/usr/bin/env python3
"""Fig 8: pristine vs -0.8 % 의 바일 4점 위치 추이.

(fig 6 은 `fig6_wcc.py` 로 분리했다 — Simphony 자체 .gnu 와 같은 내용으로
WCC 합만 그린다.)
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt

# ---------------- Fig 8 : Weyl 위치 추이 (표준 a,b,c 좌표)
RED, BLU = '#c0392b', '#1f5fd0'
pri = (0.1464927, 0.0708493)      # parent_pristine  (a,b)
m08 = (0.1294180, 0.0637718)      # -0.8%  (native (0.06377,0,0.12942) -> role a=idx2, b=idx0)
dka, dkb = m08[0] - pri[0], m08[1] - pri[1]
dk = np.hypot(dka, dkb)

fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.4, 6.4),
                             gridspec_kw=dict(width_ratios=[1.05, 1.0], wspace=0.30))

# --- (a) 4점 전체
for sx in (1, -1):
    for sy in (1, -1):
        c = RED if sx * sy > 0 else BLU
        ax.annotate('', xy=(sx * m08[0], sy * m08[1]), xytext=(sx * pri[0], sy * pri[1]),
                    arrowprops=dict(arrowstyle='-|>', color='#2b2b2b', lw=2.0,
                                    shrinkA=10, shrinkB=9), zorder=3)
        ax.plot(sx * pri[0], sy * pri[1], 'o', ms=15, color=c, mec='k', mew=1.1, zorder=5)
        ax.plot(sx * m08[0], sy * m08[1], 's', ms=13, color=c, mec='k', mew=1.1,
                alpha=0.55, zorder=4)
ax.axhline(0, color=GREY, lw=1.2, ls='--'); ax.axvline(0, color=GREY, lw=1.2, ls='--')
for s_, t_, col in ((+1, '$\\chi=+1$', RED), (-1, '$\\chi=-1$', BLU)):
    ax.text(+0.148 * s_, 0.094, t_, color=col, ha='center', fontsize=15, fontweight='bold')
    ax.text(-0.148 * s_, -0.104, t_, color=col, ha='center', fontsize=15, fontweight='bold')
ax.text(0.0, 0.122, '$k_a=0$ and $k_b=0$ are mirror lines  ($\\chi$ = 0 enforced)',
        ha='center', fontsize=12, color=GREY)
ax.set_xlim(-0.205, 0.205); ax.set_ylim(-0.135, 0.135)
ax.set_xlabel('$k_a$  (reduced)'); ax.set_ylabel('$k_b$  (reduced)')
ax.set_title('(a)  Weyl quartet in the $k_c$ = 0 plane', fontsize=15, pad=10)
ax.grid(alpha=0.18, lw=0.7)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], marker='o', ls='', ms=12, color=RED, mec='k', label='unstrained  $\\chi=+1$'),
                   Line2D([], [], marker='o', ls='', ms=12, color=BLU, mec='k', label='unstrained  $\\chi=-1$'),
                   Line2D([], [], marker='s', ls='', ms=11, color='#8a8a8a', mec='k', label='0.8 % compressive')],
          loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3,
          frameon=False, handletextpad=0.4, columnspacing=1.4)

# --- (b) 1사분면 확대
bx.annotate('', xy=m08, xytext=pri,
            arrowprops=dict(arrowstyle='-|>', color='#2b2b2b', lw=2.4,
                            shrinkA=12, shrinkB=11), zorder=3)
bx.plot(*pri, 'o', ms=18, color=RED, mec='k', mew=1.2, zorder=5)
bx.plot(*m08, 's', ms=16, color=RED, mec='k', mew=1.2, alpha=0.55, zorder=4)
bx.annotate('unstrained\n(%.5f, %.5f)' % pri, xy=pri, xytext=(-10, 16),
            textcoords='offset points', ha='right', va='bottom', fontsize=12.5)
bx.annotate('0.8 %% compressive\n(%.5f, %.5f)' % m08, xy=m08, xytext=(12, -14),
            textcoords='offset points', ha='left', va='top', fontsize=12.5)
bx.text(0.975, 0.955,
        '$|\\Delta k| = %.4f$\n$\\Delta k_a = %+.4f$\n$\\Delta k_b = %+.4f$' % (dk, dka, dkb),
        transform=bx.transAxes, ha='right', va='top', fontsize=13.5,
        bbox=dict(boxstyle='round,pad=0.45', fc='white', ec='#bbbbbb', lw=1.0, alpha=0.95))
bx.text(0.03, 0.045, 'toward $k_a=0,\\ k_b=0$\n(annihilation between $-0.8$ % and $-1.0$ %)',
        transform=bx.transAxes, ha='left', va='bottom', fontsize=12, color=GREY)
bx.set_xlim(0.1195, 0.1585); bx.set_ylim(0.0578, 0.0788)
bx.set_xlabel('$k_a$  (reduced)'); bx.set_ylabel('$k_b$  (reduced)')
bx.set_title('(b)  first quadrant, magnified', fontsize=15, pad=10)
bx.grid(alpha=0.25, lw=0.7)

fig.suptitle('Weyl motion under biaxial compression — toward the mirror lines, where the quartet annihilates',
             y=0.99, fontsize=15.5)
fig.savefig('figs/fig8_weyl_shift.png')
plt.close(fig)
print('fig8 저장')
print('이동량 (표준좌표): Delta k_a=%+.5f  Delta k_b=%+.5f  |Delta k|=%.5f' % (dka, dkb, dk))
