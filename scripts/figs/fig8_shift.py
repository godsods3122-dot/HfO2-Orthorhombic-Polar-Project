#!/usr/bin/env python3
"""Fig 8: 압축 → 인장 biaxial 에서 바일 4점 위치 추이.

⚠️ 라벨 주의 — 디렉토리 이름의 공칭값과 실제 변형률이 다르다. 역할 기준
(a = 긴 거울축, b = 짧은 거울축, c = 편극축) 실측값은 이렇다:

  m1_mirror      εa −0.687 %  εb −0.844 %  εc +0.514 %   ← 실제로는 −0.8 %
  pristine_mirror  0            0            0
  p1_mirror      εa +1.320 %  εb +1.159 %  εc −0.754 %
  parent_m1_old  εa −1.079 %  εb −0.911 %  εc +0.740 %   ← 이쪽이 진짜 −1 %

그래서 이 그림은 **−0.8 % → 0 → +1.2 %** 다. 진짜 −1 % (parent_m1_old) 에서는
polar=0 평면의 4점 궤도가 **이미 소멸했다** (201×201 스캔, 거울선 밖 국소최소
2개뿐이고 gap 7.9e-03 / 8.7e-03 — 노드보다 네 자리 크다).

좌표는 표준 (a, b, c) — 편극축 = c. pristine/m1/p1 세 구조 모두 축 규약이 같아
(role2native = [2,0,1], polar = native 1) 같은 평면에 그대로 겹쳐 그릴 수 있다.

위치는 전부 이 세션에서 다시 찾은 값이다 (기록값을 그대로 쓰지 않았다):
polar=0 평면 201x201 스캔 -> 거울선 밖 국소최소 **전부** 격자 미세화 -> 일반위치에
남는 것만 채택. 세 구조 모두 평면 위 노드는 사분면당 정확히 하나였다.
자세한 내용은 results/strain_trend/SUMMARY.md.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# 표준 (k_a, k_b), k_c = 0.  k_a = native k3, k_b = native k1.
M1 = (0.1294182, 0.0637717)      # −1 %  (native 0.0637717, 0, 0.1294182)   10.365955 THz
PRI = (0.1610282, 0.0975234)     #  0 %  (native 0.0975234, 0, 0.1610282)   10.180597 THz
P1 = (0.3182116, 0.2643381)      # +1 %  (native 0.2643381, 0, 0.3182116)    9.960553 THz
E = {'m1': 10.365955, 'pri': 10.180597, 'p1': 9.960553}

dka, dkb = P1[0] - M1[0], P1[1] - M1[1]
dk = np.hypot(dka, dkb)

C_M1, C_P1, C_PRI = '#1f5fd0', '#c0392b', '#7a7a7a'

fig, (ax, bx) = plt.subplots(1, 2, figsize=(14.8, 6.6),
                             gridspec_kw=dict(width_ratios=[1.05, 1.0], wspace=0.30))

# --- (a) 4점 전체 -------------------------------------------------------------
for sx in (1, -1):
    for sy in (1, -1):
        ax.annotate('', xy=(sx * P1[0], sy * P1[1]), xytext=(sx * M1[0], sy * M1[1]),
                    arrowprops=dict(arrowstyle='-|>', color='#2b2b2b', lw=1.8,
                                    shrinkA=9, shrinkB=9, alpha=0.75), zorder=3)
        ax.plot(sx * PRI[0], sy * PRI[1], 'o', ms=10, color='white', mec=C_PRI,
                mew=2.0, zorder=4)
        ax.plot(sx * M1[0], sy * M1[1], 'o', ms=14, color=C_M1, mec='k', mew=1.1, zorder=5)
        ax.plot(sx * P1[0], sy * P1[1], 's', ms=13, color=C_P1, mec='k', mew=1.1, zorder=5)
ax.axhline(0, color=GREY, lw=1.2, ls='--'); ax.axvline(0, color=GREY, lw=1.2, ls='--')
ax.text(0.0, -0.375, '$k_a=0$ and $k_b=0$ are mirror lines  ($\\chi$ = 0 enforced)',
        ha='center', va='top', fontsize=12, color=GREY)
ax.set_xlim(-0.43, 0.43); ax.set_ylim(-0.40, 0.36)
ax.set_xlabel('$k_a$  (reduced)'); ax.set_ylabel('$k_b$  (reduced)')
ax.set_title('(a)  Weyl quartet in the $k_c$ = 0 plane', fontsize=15, pad=14)
ax.grid(alpha=0.18, lw=0.7)
ax.set_aspect('equal')
ax.legend(handles=[
    Line2D([], [], marker='o', ls='', ms=12, color=C_M1, mec='k',
           label='compressive  $-0.84/-0.69$ %'),
    Line2D([], [], marker='o', ls='', ms=9, color='white', mec=C_PRI, mew=2.0,
           label='unstrained'),
    Line2D([], [], marker='s', ls='', ms=11, color=C_P1, mec='k',
           label='tensile  $+1.16/+1.32$ %')],
    loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=3,
    frameon=False, handletextpad=0.4, columnspacing=1.6)

# --- (b) 1사분면 확대 ---------------------------------------------------------
bx.annotate('', xy=P1, xytext=M1,
            arrowprops=dict(arrowstyle='-|>', color='#2b2b2b', lw=2.4,
                            shrinkA=13, shrinkB=12), zorder=3)
bx.plot(*M1, 'o', ms=17, color=C_M1, mec='k', mew=1.2, zorder=5)
bx.plot(*PRI, 'o', ms=12, color='white', mec=C_PRI, mew=2.4, zorder=5)
bx.plot(*P1, 's', ms=15, color=C_P1, mec='k', mew=1.2, zorder=5)
bx.annotate('$-0.84/-0.69$ %%\n(%.5f, %.5f)\n%.3f THz' % (M1 + (E['m1'],)), xy=M1,
            xytext=(16, -4), textcoords='offset points', ha='left', va='top',
            fontsize=12.5, color=C_M1)
bx.annotate('unstrained\n(%.5f, %.5f)\n%.3f THz' % (PRI + (E['pri'],)), xy=PRI,
            xytext=(16, -4), textcoords='offset points', ha='left', va='top',
            fontsize=12, color='#5a5a5a')
bx.annotate('$+1.16/+1.32$ %%\n(%.5f, %.5f)\n%.3f THz' % (P1 + (E['p1'],)), xy=P1,
            xytext=(-22, -14), textcoords='offset points', ha='right', va='top',
            fontsize=12.5, color=C_P1)
bx.text(0.975, 0.055,
        '$|\\Delta k| = %.4f$\n$\\Delta k_a = %+.4f$\n$\\Delta k_b = %+.4f$' % (dk, dka, dkb),
        transform=bx.transAxes, ha='right', va='bottom', fontsize=13.5, zorder=8,
        bbox=dict(boxstyle='round,pad=0.45', fc='white', ec='#bbbbbb', lw=1.0, alpha=0.95))
bx.set_xlim(0.100, 0.360); bx.set_ylim(0.020, 0.300)
bx.set_xlabel('$k_a$  (reduced)'); bx.set_ylabel('$k_b$  (reduced)')
bx.set_title('(b)  first quadrant, magnified', fontsize=15, pad=14)
bx.grid(alpha=0.25, lw=0.7)

fig.suptitle('Weyl position under biaxial strain  (mirror-axis strains, not the nominal labels)',
             y=0.99, fontsize=15.5)
fig.savefig('figs/fig8_weyl_shift.png')
plt.close(fig)
print('fig8 저장')
print('−1%% -> +1%%  Δk_a=%+.5f  Δk_b=%+.5f  |Δk|=%.5f' % (dka, dkb, dk))
