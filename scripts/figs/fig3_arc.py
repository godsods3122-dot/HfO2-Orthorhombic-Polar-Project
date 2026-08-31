#!/usr/bin/env python3
"""Fig 3: E = 10.0869 THz 등주파수면의 표면 상태 (left / right / bulk).

Simphony SlabArc_calc 결과 arc.dat_l / arc.dat_r / arc.dat_bulk 를 읽는다.
표면 법선은 c(편극축)이라 표면 BZ 좌표가 (k_a, k_b) 그대로다.
파일 좌표는 Cartesian 이므로 KPLANE_SLAB 로 지정한 reduced 범위에
선형 대응시켜 되돌린다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

D = sys.argv[1] if len(sys.argv) > 1 else 'work/slabarc'
KA = (-0.30, 0.30)          # KPLANE_SLAB 로 준 reduced 범위
KB = (-0.20, 0.20)
NODES = [(+0.1464927, +0.0708493, +1), (+0.1464927, -0.0708493, -1),
         (-0.1464927, +0.0708493, -1), (-0.1464927, -0.0708493, +1)]
EARC = 10.0869


def load(name):
    d = np.loadtxt(os.path.join(D, name))
    kx, ky, v = d[:, 0], d[:, 1], d[:, 2]
    # 어느 쪽이 안쪽 루프인지 데이터에서 직접 판정한다
    if abs(ky[1] - ky[0]) > abs(kx[1] - kx[0]):        # ky 가 안쪽
        n_in = int(np.argmax(np.diff(ky) < 0) + 1)
        n_out = len(v) // n_in
        V = v.reshape(n_out, n_in)                      # [kx, ky]
    else:                                               # kx 가 안쪽
        n_in = int(np.argmax(np.diff(kx) < 0) + 1)
        n_out = len(v) // n_in
        V = v.reshape(n_out, n_in).T                    # -> [kx, ky]
        n_in, n_out = n_out, n_in
    a = np.linspace(KA[0], KA[1], V.shape[0])
    b = np.linspace(KB[0], KB[1], V.shape[1])
    return a, b, V


titles = ['(a)  bottom surface', '(b)  top surface', '(c)  bulk projection']
files = ['arc.dat_l', 'arc.dat_r', 'arc.dat_bulk']
fig, axs = plt.subplots(1, 3, figsize=(17.4, 5.6), sharex=True, sharey=True)
for ax, fn, tt in zip(axs, files, titles):
    a, b, V = load(fn)
    lo, hi = np.percentile(V, [3, 99.5])
    im = ax.pcolormesh(a, b, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    for ka, kb, c in NODES:
        ax.plot(ka, kb, 'o', ms=12, mfc='none', mec='white', mew=3.0, zorder=6)
        ax.plot(ka, kb, 'o', ms=12, mfc='none', mec='k', mew=1.3, zorder=7)
        ax.text(ka, kb + 0.026, '$\\chi=%+d$' % c, color='white', ha='center',
                fontsize=13, fontweight='bold', zorder=8,
                path_effects=[pe.withStroke(linewidth=2.6, foreground='black')])
    ax.set_xlabel('$k_a$  (reduced)')
    ax.set_title(tt, fontsize=15, pad=8)
    ax.set_xlim(*KA); ax.set_ylim(*KB)
    ax.set_aspect('auto')
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label('log DOS', fontsize=12)
axs[0].set_ylabel('$k_b$  (reduced)')
fig.suptitle('Iso-frequency surface spectrum at the Weyl energy  $\\omega$ = %.4f THz  '
             '(surface $\\perp$ polar axis $c$,  $N_p$=2,  $\\eta$=0.008 THz)' % EARC,
             y=1.005, fontsize=15)
fig.savefig('figs/fig3_fermi_arc.png')
print('fig3 저장')
