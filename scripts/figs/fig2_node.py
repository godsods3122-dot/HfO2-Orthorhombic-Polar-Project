#!/usr/bin/env python3
"""Fig 2: band 17/18 하이라이트 + 실제 Weyl 노드 통과 선.

주의 1 — 안쪽 Γ 의 불연속: LO-TO 비해석항은 Γ 접근 방향에 의존하므로 Y→Γ 와
  Γ→Z 의 Γ 값이 다르다 (band17 0.0265, band18 0.0363 THz 점프). 이건 물리이지
  수치 오차가 아니다. 이어 그리면 거의 수직인 가짜 선이 생기므로 NaN 으로 끊는다.

주의 2 — Γ-X 위 교점: 격자 데이터(.dat)의 최소 gap 은 1.9e-3 THz 지만 이건 격자가
  교점을 비껴간 것이다. phonopy 로 정밀화하면 t = 0.185742 에서 gap 1.6e-11,
  E = 10.192507 THz. Γ-X 는 거울면 위 선이라 대칭이 보장하는 nodal line 점이다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY, BLUE
setup()
import matplotlib.pyplot as plt

NB, NK = 36, 120
LAB = ['Γ', 'X', 'S', 'Y', 'Γ', 'Z', 'U', 'R', 'T', 'Z']
GAMMA_BREAK = 4 * NK            # 안쪽 Γ — LO-TO 방향 의존으로 불연속

# Weyl 노드 에너지 — Γ-X 구간에 이 높이를 반투명 띠로 표시한다 (노드 자체는 경로 밖)
WEYL_E = 10.086936
# Γ-X 위에서 밴드 17/18 이 실제로 만나는 점 (phonopy 정밀화, 이 파일 docstring 참조)
XC_FRAC, XC_E = 0.185742 / 0.5, 10.192507
RED = '#c0392b'

d = np.loadtxt('figs/bulkek_parent_pristine.dat')
E = d[:, 1].reshape(NB, -1); x = d[:len(d) // NB, 0]
ticks = [x[i * NK] for i in range(len(LAB) - 1)] + [x[-1]]

# 불연속 지점에 NaN 을 끼워 선을 끊는다
xp = np.insert(x, GAMMA_BREAK, np.nan)
Ep = np.insert(E, GAMMA_BREAK, np.nan, axis=1)

c = np.loadtxt('figs/bulkek_nodecut.dat')
Ec = c[:, 1].reshape(NB, -1); xc = c[:len(c) // NB, 0]
n = Ec.shape[1] // 2
xa, Ea = xc[:n], Ec[:, :n]
xb, Eb = xc[n:] - xc[n], Ec[:, n:]

fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.6),
                         gridspec_kw=dict(width_ratios=[2.15, 1, 1], wspace=0.28))

ax = axes[0]
for b in range(NB):
    if b not in (16, 17):
        ax.plot(xp, Ep[b], color='#c8c8c8', lw=1.0)
# Γ-X 구간에만 Weyl 노드 에너지를 반투명 굵은 선으로
ax.plot([ticks[0], ticks[1]], [WEYL_E, WEYL_E], color=RED, lw=5.5, alpha=0.32,
        solid_capstyle='butt', zorder=1)
ax.plot(xp, Ep[16], color=BLUE, lw=2.5, label='band 17', zorder=3)
ax.plot(xp, Ep[17], color='#67b0ff', lw=2.5, label='band 18', zorder=3)
for t in ticks[1:-1]:
    ax.axvline(t, color=GREY, lw=0.8, alpha=0.5)
# 두 밴드가 만나는 점
xcross = ticks[0] + XC_FRAC * (ticks[1] - ticks[0])
ax.plot(xcross, XC_E, 'o', ms=10, mfc='none', mec=RED, mew=2.4, zorder=4)

ax.set_xticks(ticks)
ax.set_xticklabels(LAB, fontsize=15, fontweight='bold')   # 고대칭 경로 볼드
ax.set_xlim(x[0], x[-1])
ax.set_ylim(9.1, 11.95); ax.set_ylabel('Frequency (THz)')
from matplotlib.lines import Line2D
h, l = ax.get_legend_handles_labels()
h += [Line2D([], [], color=RED, lw=5.5, alpha=0.32),
      Line2D([], [], ls='none', marker='o', ms=9, mfc='none', mec=RED, mew=2.2)]
l += ['Weyl node  %.3f THz' % WEYL_E,
      'bands 17/18 meet on $\\Gamma$–X  %.3f THz' % XC_E]
ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, -0.155), ncol=2,
          frameon=False, fontsize=11.5, handlelength=2.4,
          columnspacing=2.4, labelspacing=0.7)
ax.set_title('(a)  bands 17 / 18 on the standard path', loc='left')

for ax, xx, EE, ttl, sub in ((axes[1], xa, Ea, '(b)  cut along $k_a$', '$k_b$=0.0708, $k_c$=0'),
                             (axes[2], xb, Eb, '(c)  cut along $k_b$', '$k_a$=0.1465, $k_c$=0')):
    for b in range(NB):
        if b not in (16, 17):
            ax.plot(xx, EE[b], color='#c8c8c8', lw=1.0)
    # (a) 와 같은 높이의 Weyl 노드 에너지 띠 — 여기서는 교점이 곧 노드다
    ax.axhline(WEYL_E, color=RED, lw=5.5, alpha=0.32, zorder=1)
    ax.plot(xx, EE[16], color=BLUE, lw=2.6, zorder=3)
    ax.plot(xx, EE[17], color='#67b0ff', lw=2.6, zorder=3)
    g = EE[17] - EE[16]; i = int(np.argmin(g))
    ax.plot(xx[i], EE[16][i], 'o', ms=9, mfc='none', mec=RED, mew=2.2, zorder=4)
    ax.set_xticks([xx[0], xx[i], xx[-1]]); ax.set_xticklabels(['−0.035', '0', '+0.035'])
    ax.set_xlabel('Δk (reduced)')
    ax.set_ylim(EE[16].min() - 0.09, EE[17].max() + 0.09)
    ax.set_title(ttl, loc='left'); ax.text(0.5, 0.03, sub, transform=ax.transAxes,
                                           ha='center', fontsize=11, color='#555')
axes[1].set_ylabel('Frequency (THz)')
fig.text(0.655, 0.955, 'Weyl node  $k$=(0.1465, 0.0708, 0),  10.087 THz = 41.72 meV,  $\\chi=+1$',
         ha='center', fontsize=12.5, color=RED)
fig.savefig('figs/fig2_bands17_18_node.png')
print('fig2 저장. 노드 선 위 최소 gap = %.2e THz' % (Ec[17] - Ec[16]).min())
print('Γ-X 교점 E = %.6f THz, x 위치 = %.4f' % (XC_E, xcross))
