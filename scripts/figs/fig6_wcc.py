#!/usr/bin/env python3
"""Fig 6: 네 바일 노드의 구면 Wilson loop (WCC).

Simphony `WeylChirality_calc` 가 뽑아 주는 `wanniercenter3D_Weyl_{1..4}.gnu` 와
**같은 내용**이다.  그 스크립트들은 `wanniercenter3D_Weyl.dat` 의
열 2, 3, 4, 5 (= 네 노드 각각의 WCC 합) 만 South→North 로 찍는다.
합의 감김수가 곧 chirality 이므로, 밴드별 WCC 산점도를 같이 그리면
오히려 가려서 안 보인다 — 그래서 합만 그린다.

이 컨테이너에 gnuplot 이 없어 같은 열·같은 축 범위로 matplotlib 에 옮겼다.
(`work/ppchi/wanniercenter3D_Weyl_1.gnu` 참조)
"""
import sys, os, re, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else 'work/ppchi'
d = np.loadtxt(os.path.join(D, 'wanniercenter3D_Weyl.dat'))
t = d[:, 0]

# Simphony 가 만든 .gnu 에서 노드 좌표와 사용 열 번호를 그대로 읽어 온다
info = []
for g in sorted(glob.glob(os.path.join(D, 'wanniercenter3D_Weyl_*.gnu'))):
    s = open(g).read()
    ttl = re.search(r'set title "Weyl \(([^)]*)\)', s).group(1)
    col = int(re.search(r'u 1:\s*(\d+)', s).group(1))
    info.append((ttl.strip(), col))


def winding(y):
    """WCC 합의 감김수 = 언랩한 증분의 총합 (반올림)."""
    u = np.unwrap(np.asarray(y) * 2 * np.pi) / (2 * np.pi)
    return int(round(u[-1] - u[0]))


RED, BLU = '#c0392b', '#1f5fd0'
fig, axs = plt.subplots(1, 4, figsize=(17.6, 4.9), sharey=True,
                        gridspec_kw=dict(wspace=0.16))
for ax, (ttl, col) in zip(axs, info):
    y = d[:, col - 1]
    w = winding(y)
    c = RED if w > 0 else BLU
    ax.plot(t, y, 'o', ms=2.6, color=c, rasterized=True)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(['S', '', 'equator', '', 'N'])
    ax.set_xlabel('$k$  on the sphere')
    ka, kb, kc = [float(x) for x in ttl.split(',')]
    ax.set_title('$k$ = (%+.5f, %+.5f, %.0f)\n$\\chi = %+d$' % (ka, kb, kc, w),
                 fontsize=13, pad=8, color=c, fontweight='bold')
    ax.grid(alpha=0.22, lw=0.7)
axs[0].set_ylabel('WCC')
fig.suptitle('Wilson loops on spheres around the four Weyl nodes '
             '— the winding of the WCC sum is the chirality\n'
             '$r_0$ = 0.0015 $\\AA^{-1}$,  $N_{k1}$ = 240,  $N_{k2}$ = 801,  bands 1–17   '
             '(same quantity as Simphony\'s own wanniercenter3D_Weyl_*.gnu)',
             y=1.13, fontsize=14)
fig.savefig('figs/fig6_wcc_4weyl.png')
print('fig6 저장.  chirality =', [winding(d[:, c - 1]) for _, c in info])
