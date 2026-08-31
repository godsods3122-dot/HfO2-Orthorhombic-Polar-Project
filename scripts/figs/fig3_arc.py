#!/usr/bin/env python3
"""Fig 3: ω = 10.0869 THz (바일 주파수) 에서의 표면 스펙트럼 한 장.

`dos_l − dos_r` (아래 표면 − 위 표면).  벌크 기여가 두 표면에 똑같이 들어가
정확히 상쇄되므로 표면 성분만 남는다.  Simphony `SlabArc_calc` 의
`arc.dat_l` / `arc.dat_r` 를 읽는다.  표면 법선 ∥ c (편극축) 이라 표면 BZ
좌표가 (k_a, k_b) 그대로다.

⚠️ 이 그림에 보이는 밝은 능선은 **페르미 아크가 아니다.**  능선을 추적한 결과:

  k_b        능선 k_a    배경 대비
  0.0000     0.0975      1.36
  0.0708     0.0950      1.30      <- 노드의 k_b.  아무 일도 일어나지 않는다
  0.1000     0.0875      0.99
  0.1200     (다른 능선)  0.68

노드의 k_b 를 그냥 통과해 |k_b| ≈ 0.11 까지 이어지다 사라진다.  아크라면
사영점에서 끝나는 **열린** 곡선이어야 하는데 그렇지 않다.  게다가 이 능선이
있는 k_a ≈ 0.0975 는 band 17 과 18 이 **둘 다** 사영된 영역이다.  아크가
존재할 수 있는 곳은 둘 중 하나만 사영된 영역인데, 거기서는(k_b=0 기준
k_a 0.1175~0.1700) dos_l−dos_r 최대가 0.563 으로 배경 중앙값 0.444 와
구분되지 않는다.

즉 이 계산에서 **식별 가능한 페르미 아크는 없다.**  근본 원인은 노드
주파수에 벌크 간극이 없다는 것이다 — 사영 벌크 연속체가 표면 BZ 창 전체를
덮는다 (phonopy 확인, 커버리지 1.000).  36개 가지가 이 대역에 빽빽해서
아크가 존재하더라도 연속체에 묻힌 공명이다.

자세한 내용은 `results/parent_pristine/SUMMARY.md` 17절.
"""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

DIRS = sys.argv[1:] or sorted(glob.glob('work/arc_loto_*'))
XL, YL = (-0.30, 0.30), (-0.20, 0.20)
NODES = [(+0.1464927, +0.0708493, +1), (+0.1464927, -0.0708493, -1),
         (-0.1464927, +0.0708493, -1), (-0.1464927, -0.0708493, +1)]
EARC = 10.0869
STK = [pe.withStroke(linewidth=3.4, foreground='black')]


def load_one(d, name):
    a = np.loadtxt(os.path.join(d, name))
    kx, ky, v = a[:, 0], a[:, 1], a[:, 2]
    n_in = int(np.argmax(np.diff(ky) < 0) + 1)          # ky 가 안쪽 루프
    n_out = len(v) // n_in
    s = open(os.path.join(d, 'PN.out')).read()
    st = [float(x) for x in re.search(r'K2D_start:\s*(\S+)\s+(\S+)', s).groups()]
    v1 = [float(x) for x in re.search(r'The first vector:\s*(\S+)\s+(\S+)', s).groups()]
    v2 = [float(x) for x in re.search(r'The second vector:\s*(\S+)\s+(\S+)', s).groups()]
    return (st[0] + v1[0] * np.linspace(0, 1, n_out),
            st[1] + v2[1] * np.linspace(0, 1, n_in), v.reshape(n_out, n_in))


def stitch(name):
    """k_a 로 쪼개 병렬로 돌린 조각을 붙인다 (경계 중복 열은 버린다)."""
    P = [load_one(d, name) for d in DIRS]
    P.sort(key=lambda p: p[0][0])
    ka, kb, V = P[0][0], P[0][1], P[0][2]
    for k2, _, V2 in P[1:]:
        if abs(k2[0] - ka[-1]) < 1e-9:
            k2, V2 = k2[1:], V2[1:]
        ka = np.concatenate([ka, k2]); V = np.concatenate([V, V2], axis=0)
    return ka, kb, V


a, b, L = stitch('arc.dat_l')
_, _, R = stitch('arc.dat_r')

# 표면 BZ 거울 대칭을 명시적으로 회복시킨다.  이 표면(법선 ∥ c)은 z 이동이 없는
# a-glide(⊥b)를 그대로 가지므로 dos(k_a,k_b) = dos(k_a,-k_b) 여야 한다.
# 시간역전에서 오는 반전 대칭은 계산값이 이미 정확히(0.000e+00) 만족하는데
# 거울은 log DOS 중앙값 0.024 / 최대 0.66 (범위 7.8) 어긋난다.  벌크 phonopy 는
# 같은 거울을 1e-7 THz 로 지키므로 입력이 아니라 ham_qlayer2qlayer_LOTO 쪽
# 문제로 보인다 (SUMMARY 17.6).  원인을 못 짚었으므로 알려진 대칭으로 평균만 낸다.
L, R = (0.5 * (V + V[:, ::-1]) for V in (L, R))
V = L - R

fig, ax = plt.subplots(figsize=(9.2, 7.4))
m = (a >= XL[0]) & (a <= XL[1]); n = (b >= YL[0]) & (b <= YL[1])
lo, hi = np.percentile(V[np.ix_(m, n)], [35, 99.6])
im = ax.pcolormesh(a, b, V.T, cmap='jet', vmin=lo, vmax=hi,
                   shading='gouraud', rasterized=True)
for ka, kb, c in NODES:
    ax.plot(ka, kb, 'o', ms=14, mfc='none', mec='white', mew=3.4, zorder=6)
    ax.plot(ka, kb, 'o', ms=14, mfc='none', mec='k', mew=1.4, zorder=7)
    ax.annotate('$\\chi=%+d$' % c, xy=(ka, kb),
                xytext=(0, 26 if kb > 0 else -26), textcoords='offset points',
                color='white', ha='center', va='center', fontsize=14,
                fontweight='bold', zorder=8, path_effects=STK)
ax.set_xlim(*XL); ax.set_ylim(*YL)
ax.set_xlabel('$k_a$  (reduced)')
ax.set_ylabel('$k_b$  (reduced)')
cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
cb.set_label('$\\log(\\rho_{\\rm bottom}/\\rho_{\\rm top})$', fontsize=13)
ax.set_title('Surface spectral weight at the Weyl frequency\n'
             '$\\omega$ = %.4f THz,  surface $\\perp$ polar axis $c$' % EARC,
             fontsize=15, pad=12)
fig.savefig('figs/fig3_surface_spectrum.png')
print('fig3 저장.  격자 %d x %d,  간격 %.5f x %.5f'
      % (len(a), len(b), a[1] - a[0], b[1] - b[0]))
