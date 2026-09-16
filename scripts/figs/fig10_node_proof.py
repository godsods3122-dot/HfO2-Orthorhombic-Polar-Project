#!/usr/bin/env python3
"""노드가 진짜인지 눈으로 확인하는 그림: 고립된 점 영점 + 선형 V자.

(a) 노드를 지나는 평면의 band17-18 gap 지도 (로그, 단일 색상 순차 램프)
(b) 노드를 지나는 세 방향 1D 절단 — 선형이면 V자, 2차 접촉이면 포물선

사용: python3 fig10_node_proof.py SRC K1 K2 K3 OUT.png [TITLE]
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
sys.path.insert(0, __file__.rsplit('/', 2)[0])
from style import setup, GREY
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from weyl_scan import get_ph, gap_at

SRC = sys.argv[1]
K0 = np.array([float(x) for x in sys.argv[2:5]])
OUT = sys.argv[5]
TITLE = sys.argv[6] if len(sys.argv) > 6 else ''
BAND = 17
setup()

ph = get_ph(SRC)
g = lambda Q: gap_at(ph, [list(q) for q in Q], BAND)

# --- (a) 콘이 가장 잘 열리는 두 축을 고른다: 속도행렬 특이값이 큰 쪽 -------------
W = 0.02
N = 121
ax = np.linspace(-W, W, N)
A, B = np.meshgrid(ax, ax, indexing='ij')
I, J = 0, 2                                 # k1, k3 평면
Q = np.tile(K0, (N * N, 1))
Q[:, I] += A.ravel(); Q[:, J] += B.ravel()
G = g(Q).reshape(N, N)

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4))
axp = axes[0]
vmin = max(G.min(), 1e-6)
im = axp.pcolormesh(K0[I] + ax, K0[J] + ax, G.T, cmap='Blues',
                    norm=LogNorm(vmin, G.max()), shading='auto', rasterized=True)
axp.plot(K0[I], K0[J], 'o', mfc='none', mec='#b32424', mew=2.4, ms=15)
axp.annotate('node', (K0[I], K0[J]), textcoords='offset points', xytext=(16, 12),
             color='#b32424', fontsize=13, fontweight='bold')
axp.set_xlabel('$k_%d$ (reduced)' % (I + 1))
axp.set_ylabel('$k_%d$ (reduced)' % (J + 1))
axp.set_title('(a)  band 17–18 gap, $k_%d$ fixed' % (3 - I - J + 1), loc='left')
axp.set_aspect('equal')
cb = fig.colorbar(im, ax=axp); cb.set_label('gap (THz)')

# --- (b) 세 방향 1D 절단 -------------------------------------------------------
axc = axes[1]
t = np.concatenate([-np.logspace(-4.5, np.log10(W), 60)[::-1], [0],
                    np.logspace(-4.5, np.log10(W), 60)])
cols = ['#1f5fd0', '#d17a00', '#1a7f4f']
for i, c in zip(range(3), cols):
    Q = np.tile(K0, (len(t), 1)); Q[:, i] += t
    axc.plot(t, g(Q), '-', lw=2.0, color=c, label='along $k_%d$' % (i + 1))
axc.axvline(0, color=GREY, lw=1.0, ls='--', alpha=0.7)
axc.set_yscale('log')
axc.set_xlabel('displacement from node (reduced)')
axc.set_ylabel('band 17–18 gap (THz)')
axc.set_title('(b)  1D cuts through the node — a linear crossing gives a V', loc='left')
axc.legend(frameon=False, loc='lower right')

if TITLE:
    fig.suptitle(TITLE, fontsize=15)
fig.tight_layout()
fig.savefig(OUT)
print('wrote', OUT, ' 평면 최소 gap %.3e THz' % G.min())
