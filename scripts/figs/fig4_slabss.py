#!/usr/bin/env python3
"""Fig 4: 페르미 아크를 가로지르는 경로 위의 표면 상태 스펙트럼.

경로는 `k_b = 0` 고정, `k_a : 0 → 0.30`.  fig 3 (d) 에서 아크는 chirality 가
반대인 두 사영점 `(0.14649, ±0.07085)` 을 잇고 `k_b = 0` 에서 `k_a ≈ 0.110`
까지 안쪽으로 휘어 있다.  그래서 이 직선은 **아크를 반드시 한 번 관통한다.**

노드 주파수에는 벌크 간극이 없으므로(사영 연속체가 표면 BZ 전체를 덮는다)
raw DOS 로는 벌크에 묻힌다.  (c) 의 표면 전용 값에서 아크 가지가 ω_Weyl 을
가로지르는 게 선명하게 보인다.

표면 전용 = ρ_surf − ρ_bulk·N_top/N_dim.  Simphony 의 `dos_l_only` 열은 이
궤도 개수 보정이 없어 전부 eps9 로 잘려 나온다
(`patch/apply_surfdos_only_norm_fix.py`).
"""
import sys, os, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

D = sys.argv[1] if len(sys.argv) > 1 else 'work/cut_zero'
KA0, KA1 = 0.0, 0.30
W0 = 10.086936
KARC = 0.110              # 아크가 이 선을 지나는 k_a (데이터에서 읽은 값)
STK = [pe.withStroke(linewidth=3.2, foreground='black')]


def load(name):
    d = np.loadtxt(os.path.join(D, name))
    k, e, v = d[:, 0], d[:, 1], d[:, 2]
    no = int(np.argmax(np.diff(e) < 0) + 1)          # omega 가 안쪽 루프
    nk = len(k) // no
    kk = k.reshape(nk, no)[:, 0]
    ka = KA0 + (KA1 - KA0) * (kk - kk[0]) / (kk[-1] - kk[0])
    return ka, e.reshape(nk, no)[0], v.reshape(nk, no)


ka, om, L = load('dos.dat_l')
_, _, R = load('dos.dat_r')
_, _, BK = load('dos.dat_bulk')

_s = open(os.path.join(D, 'PN.out')).read()
RATIO = (int(re.search(r'NtopOrbitals\s+(\d+)', _s).group(1))
         / int(re.search(r'(?i)\bndim:?\s+(\d+)', _s).group(1)))
ONLY = np.log10(np.maximum(np.exp(L) - np.exp(BK) * RATIO, 0.0) + 1e-3)

panels = [(L, '(a)  bottom surface', 'log DOS'),
          (R, '(b)  top surface', 'log DOS'),
          (ONLY, '(c)  bottom surface only  —  the Fermi arc branch',
           '$\\log_{10}$ surface DOS')]
fig, axs = plt.subplots(1, 3, figsize=(18.0, 6.0), sharey=True,
                        gridspec_kw=dict(wspace=0.16))
for ax, (V, tt, cl) in zip(axs, panels):
    lo, hi = np.percentile(V, [3, 99.4])
    im = ax.pcolormesh(ka, om, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    ax.axhline(W0, color='white', lw=1.7, ls='--', alpha=0.9, path_effects=STK)
    ax.plot(KARC, W0, 'o', ms=14, mfc='none', mec='white', mew=3.2, zorder=7)
    ax.plot(KARC, W0, 'o', ms=14, mfc='none', mec='k', mew=1.3, zorder=8)
    ax.set_xlim(KA0, KA1)
    ax.set_xlabel('$k_a$  (reduced)   at  $k_b$ = 0')
    ax.set_title(tt, fontsize=14.5, pad=8)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(cl, fontsize=11.5)
axs[0].set_ylabel('Frequency (THz)')
axs[0].text(0.025, 0.975, '$\\omega_{\\rm Weyl}$ = %.4f THz' % W0, color='white',
            transform=axs[0].transAxes, va='top', fontsize=12.5, path_effects=STK)
axs[2].annotate('arc crosses $\\omega_{\\rm Weyl}$\nat $k_a$ = %.3f' % KARC,
                xy=(KARC, W0), xytext=(0.215, 10.19), ha='center',
                color='white', fontsize=13, fontweight='bold', zorder=9,
                path_effects=STK,
                arrowprops=dict(arrowstyle='-|>', color='white', lw=2.1,
                                shrinkB=9, path_effects=STK))
fig.suptitle('Surface spectrum on a line that must cross the Fermi arc  '
             '($k_b$ = 0,  surface $\\perp$ $c$,  LO-TO on,  $N_p$ = 2,  $\\eta$ = 0.004 THz)\n'
             'the arc connects $(0.14649, +0.07085)$ and $(0.14649, -0.07085)$ '
             'and bows inward to $k_a \\approx 0.110$ at $k_b$ = 0',
             y=1.035, fontsize=15)
fig.savefig('figs/fig4_slab_arc.png')
print('fig4 저장.  N_top/N_dim = %.3f' % RATIO)
