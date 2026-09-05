#!/usr/bin/env python3
"""Fig 5: Weyl 콘 3D.

노드 주변에서 band 17/18 을 (k_a,k_b) 격자에 올린다.  이 노드는 k_a 방향
tilt(2.44 THz/rlu)가 콘 기울기(0.97 THz/rlu)보다 커서 type-II 다 —
k_a 를 따라가면 두 가지가 모두 위로 올라간다.  등주파수면이 점이 아니라
서로 맞닿은 두 개의 포켓이 되는 게 type-II 의 지문이고, 그게 페르미 아크
그림(fig 3)의 배경 무늬로 그대로 나타난다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt
from weyl_scan import get_ph

W = np.array([0.1464927, 0.0708493, 0.0])
R, N = 0.0085, 101
ph = get_ph('source/parent_pristine')
a = np.linspace(-R, R, N)
ph.run_qpoints([[W[0] + x, W[1] + y, 0.0] for x in a for y in a])
f = ph.qpoints.frequencies
E17 = f[:, 16].reshape(N, N)
E18 = f[:, 17].reshape(N, N)
A, B = np.meshgrid(a, a, indexing='ij')
i = np.unravel_index(np.argmin(E18 - E17), E17.shape)
w0 = E17[i]

BLU, ORG = '#2b6cb0', '#dd6b20'
fig = plt.figure(figsize=(15.2, 7.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.26,
                      left=0.015, right=0.975, bottom=0.17, top=0.92)

# ---------------- (a) 3D 콘
ax = fig.add_subplot(gs[0], projection='3d')
ax.plot_surface(A, B, E18, color=ORG, alpha=0.62, linewidth=0, antialiased=True,
                rstride=2, cstride=2, shade=True)
ax.plot_surface(A, B, E17, color=BLU, alpha=0.62, linewidth=0, antialiased=True,
                rstride=2, cstride=2, shade=True)
# 노드를 지나는 두 주축 절단선을 굵게 → 교차점이 눈에 보인다
j = i[0]
for E, c in ((E17, BLU), (E18, ORG)):
    ax.plot(a, np.full(N, 0.0), E[:, j], color=c, lw=3.4, zorder=9)
    ax.plot(np.full(N, 0.0), a, E[j, :], color=c, lw=3.4, ls='--', zorder=9)
ax.scatter([0], [0], [w0], s=95, c='k', depthshade=False, zorder=12)
ax.text2D(0.60, 0.80, 'Weyl,  $\\chi=+1$\n10.0869 THz', transform=ax.transAxes,
          fontsize=13, ha='left', va='center',
          bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#cccccc', alpha=0.92))
ax.annotate('', xy=(0.505, 0.505), xytext=(0.60, 0.775), xycoords='axes fraction',
            textcoords='axes fraction',
            arrowprops=dict(arrowstyle='-', color='#444', lw=1.4))
ax.set_xlabel('\n$\\Delta k_a$', linespacing=2.2)
ax.set_ylabel('\n$\\Delta k_b$', linespacing=2.2)
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('Frequency (THz)', labelpad=32, rotation=90)
ax.set_title('(a)  Weyl cone,  $k$ = (0.14649, 0.07085, 0)', fontsize=15, pad=-2)
ax.view_init(elev=18, azim=-61)
ax.set_box_aspect((1, 1, 0.78), zoom=1.06)
ax.set_zlim(E17.min(), E18.max())
ax.tick_params(labelsize=10.5, pad=1)
ax.tick_params(axis='z', pad=9)
from matplotlib.ticker import MaxNLocator
ax.xaxis.set_major_locator(MaxNLocator(5)); ax.yaxis.set_major_locator(MaxNLocator(5))
ax.zaxis.set_major_locator(MaxNLocator(5))
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=BLU, lw=4, label='band 17'),
                   Line2D([], [], color=ORG, lw=4, label='band 18'),
                   Line2D([], [], color='k', lw=2.4, ls='-', label='cut along $\\Delta k_a$'),
                   Line2D([], [], color='k', lw=2.4, ls='--', label='cut along $\\Delta k_b$')],
          loc='upper center', bbox_to_anchor=(0.5, -0.02), ncol=2,
          frameon=False, fontsize=12)

# ---------------- (b) type-II 증거: 두 주축 절단
bx = fig.add_subplot(gs[1])
bx.plot(a, E17[:, j], color=BLU, lw=3.0, label='band 17,  along $\\Delta k_a$')
bx.plot(a, E18[:, j], color=ORG, lw=3.0, label='band 18,  along $\\Delta k_a$')
bx.plot(a, E17[j, :], color=BLU, lw=2.4, ls='--', label='band 17,  along $\\Delta k_b$')
bx.plot(a, E18[j, :], color=ORG, lw=2.4, ls='--', label='band 18,  along $\\Delta k_b$')
bx.axhline(w0, color=GREY, lw=1.1, ls=':')
bx.plot([0], [w0], 'o', ms=11, color='k', zorder=6)
bx.set_xlim(-R, R)
bx.set_xlabel('$\\Delta k$  (reduced)')
bx.set_ylabel('Frequency (THz)')
bx.set_title('(b)  cuts through the node', fontsize=15, pad=10)
bx.grid(alpha=0.22, lw=0.7)
bx.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2,
          fontsize=11.5, frameon=False)
bx.text(0.97, 0.035,
        'along $\\Delta k_a$ both branches rise\n'
        '$\\Rightarrow$ tilt 2.44 > cone slope 0.97 THz/rlu\n'
        r'$\Rightarrow$ ' + 'type-II Weyl point',
        transform=bx.transAxes, ha='right', va='bottom', fontsize=12.5,
        bbox=dict(boxstyle='round,pad=0.45', fc='#fff8e6', ec='#d8b25a', lw=1.2))

fig.savefig('figs/fig5_weyl_cone.png', bbox_inches=None)
print('fig5 저장.  gap = %.3e THz,  w0 = %.6f THz' % ((E18 - E17)[i], w0))
