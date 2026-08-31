#!/usr/bin/env python3
"""Fig 4: 바일 사영점 하나를 감싸는 닫힌 경로 위의 표면 상태 스펙트럼.

chirality 가 ±1 이면 그 사영점을 감싸는 임의의 닫힌 고리를 페르미 아크가
반드시 홀수 번 가로지른다.  그래서 이 경로는 "아크가 반드시 있는" 경로다.
Simphony SlabSS_calc 의 dos.dat_l / dos.dat_r / dos.dat_bulk 를 읽는다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt

D = sys.argv[1] if len(sys.argv) > 1 else 'work/slabloop'
EARC = 10.086936
CORNERS = ['A', 'B', 'C', 'D', 'A']


def load(name):
    d = np.loadtxt(os.path.join(D, name))
    k = d[:, 0]; e = d[:, 1]; v = d[:, 2]
    no = int(np.argmax(e[1:] < e[:-1]) + 1)          # omega 가 안쪽 루프
    nk = len(k) // no
    return k.reshape(nk, no)[:, 0], e.reshape(nk, no)[0], v.reshape(nk, no)


titles = ['(a)  bottom surface', '(b)  top surface', '(c)  bulk projection']
files = ['dos.dat_l', 'dos.dat_r', 'dos.dat_bulk']
fig, axs = plt.subplots(1, 3, figsize=(17.4, 5.8), sharey=True)
for ax, fn, tt in zip(axs, files, titles):
    k, e, V = load(fn)
    lo, hi = np.percentile(V, [3, 99.5])
    im = ax.pcolormesh(k, e, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    ax.axhline(EARC, color='w', lw=1.6, ls='--', alpha=0.85)
    n = len(CORNERS) - 1
    for j in range(1, n):
        ax.axvline(k[0] + (k[-1] - k[0]) * j / n, color='w', lw=1.0, alpha=0.55)
    ax.set_xticks([k[0] + (k[-1] - k[0]) * j / n for j in range(n + 1)])
    ax.set_xticklabels(CORNERS)
    ax.set_xlim(k[0], k[-1])
    ax.set_title(tt, fontsize=15, pad=8)
    ax.set_xlabel('closed loop around the $\\chi=+1$ projection')
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label('log DOS', fontsize=12)
axs[0].set_ylabel('Frequency (THz)')
axs[0].text(0.02, 0.965, '$\\omega_{\\rm Weyl}$ = %.4f THz' % EARC, color='w',
            transform=axs[0].transAxes, va='top', fontsize=12.5)
fig.suptitle('Surface spectrum on a closed loop enclosing one Weyl projection — '
             'the Fermi arc must cross it an odd number of times',
             y=1.005, fontsize=15)
fig.savefig('figs/fig4_slab_arc_loop.png')
print('fig4 저장')
