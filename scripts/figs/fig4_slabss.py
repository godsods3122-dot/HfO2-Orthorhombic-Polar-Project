#!/usr/bin/env python3
"""Fig 4: k_b = 0 직선 위의 표면 상태 스펙트럼 (보조 그림).

경로는 `k_b = 0` 고정, `k_a : 0 → 0.30`.

⚠️ 이 그림의 표면 가지를 한때 페르미 아크로 봤으나 **아니다.**  능선을 2D
지도에서 추적하니 노드의 `k_b = ±0.07085` 를 그냥 통과해 `|k_b| ≈ 0.11`
까지 이어지고, 그 능선이 있는 `k_a ≈ 0.10` 은 band 17 과 18 이 둘 다 사영된
영역이다.  아크는 사영점에서 끝나는 열린 곡선이어야 하고 한 밴드만 사영된
영역에 있어야 한다.  자세한 내용은 `results/parent_pristine/SUMMARY.md` 17.5.

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
KBRANCH = 0.110           # 표면 가지가 이 선에서 omega_Weyl 을 지나는 k_a
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
          (ONLY, '(c)  bottom surface only  —  surface branch',
           '$\\log_{10}$ surface DOS')]
fig, axs = plt.subplots(1, 3, figsize=(18.0, 6.0), sharey=True,
                        gridspec_kw=dict(wspace=0.16))
for ax, (V, tt, cl) in zip(axs, panels):
    lo, hi = np.percentile(V, [3, 99.4])
    im = ax.pcolormesh(ka, om, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    ax.axhline(W0, color='white', lw=1.7, ls='--', alpha=0.9, path_effects=STK)
    ax.plot(KBRANCH, W0, 'o', ms=14, mfc='none', mec='white', mew=3.2, zorder=7)
    ax.plot(KBRANCH, W0, 'o', ms=14, mfc='none', mec='k', mew=1.3, zorder=8)
    ax.set_xlim(KA0, KA1)
    ax.set_xlabel('$k_a$  (reduced)   at  $k_b$ = 0')
    ax.set_title(tt, fontsize=14.5, pad=8)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(cl, fontsize=11.5)
axs[0].set_ylabel('Frequency (THz)')
axs[0].text(0.025, 0.975, '$\\omega_{\\rm Weyl}$ = %.4f THz' % W0, color='white',
            transform=axs[0].transAxes, va='top', fontsize=12.5, path_effects=STK)
axs[2].annotate('surface branch crosses\n$\\omega_{\\rm Weyl}$ at $k_a$ = %.3f' % KBRANCH,
                xy=(KBRANCH, W0), xytext=(0.215, 10.19), ha='center',
                color='white', fontsize=13, fontweight='bold', zorder=9,
                path_effects=STK,
                arrowprops=dict(arrowstyle='-|>', color='white', lw=2.1,
                                shrinkB=9, path_effects=STK))
fig.suptitle('Surface spectrum along $k_b$ = 0  '
             '(surface $\\perp$ $c$,  LO-TO on,  $N_p$ = 2,  $\\eta$ = 0.004 THz)\n'
             'the branch crossing $\\omega_{\\rm Weyl}$ is a surface resonance inside the bulk '
             'continuum, not a Fermi arc',
             y=1.035, fontsize=15)
fig.savefig('figs/fig4_slab_arc.png')
print('fig4 저장.  N_top/N_dim = %.3f' % RATIO)
