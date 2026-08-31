#!/usr/bin/env python3
"""Fig 4: 두 바일 사영점을 관통하는 직선 위의 표면 상태 스펙트럼.

경로는 `k_a = 0.14649` 고정, `k_b : -0.25 → +0.25`.  이 선은 chirality 가
반대인 두 사영점 `(0.14649, ±0.07085)` 을 정확히 지난다.  아크는 그 두 점에서
끝나야 하므로, 이 선을 따라가면 표면 상태가 사영점에서 벌크 연속체로
빨려 들어가며 사라지는 게 보인다.

노드 주파수에서는 사영 벌크 연속체가 표면 BZ 전체를 덮으므로(간극이 없다)
아크는 고립 상태가 아니라 공명이다.  raw DOS 로는 벌크에 묻히고,
위/아래 표면 차 `dos_l − dos_r` 에서 깨끗하게 드러난다.
"""
import sys, os, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

D = sys.argv[1] if len(sys.argv) > 1 else 'work/slabcut'
KA = 0.1464927
KBN = 0.0708493
W0 = 10.086936
KB0, KB1 = -0.25, 0.25


def load(name):
    d = np.loadtxt(os.path.join(D, name))
    k, e, v = d[:, 0], d[:, 1], d[:, 2]
    no = int(np.argmax(np.diff(e) < 0) + 1)            # omega 가 안쪽 루프
    nk = len(k) // no
    kk = k.reshape(nk, no)[:, 0]
    kb = KB0 + (KB1 - KB0) * (kk - kk[0]) / (kk[-1] - kk[0])
    return kb, e.reshape(nk, no)[0], v.reshape(nk, no)


kb, om, L = load('dos.dat_l')
_, _, R = load('dos.dat_r')
_, _, BK = load('dos.dat_bulk')

# 표면 전용 DOS.  dos_l 은 표면 단위포 궤도(NtopOrbitals)만, dos_bulk 은 주층
# 전체(Ndim = Np*Num_wann)를 합하므로 궤도 개수 비로 맞춰야 한다.  Simphony 의
# dos_l_only 열은 이 보정이 없어 전부 eps9 로 잘려 나온다
# (patch/apply_surfdos_only_norm_fix.py).
_s = open(os.path.join(D, 'PN.out')).read()
RATIO = (int(re.search(r'NtopOrbitals\s+(\d+)', _s).group(1))
         / int(re.search(r'(?i)\bndim:?\s+(\d+)', _s).group(1)))
ONLY = np.log10(np.maximum(np.exp(L) - np.exp(BK) * RATIO, 0.0) + 1e-3)

STK = [pe.withStroke(linewidth=3.2, foreground='black')]
panels = [(L, '(a)  bottom surface', 'log DOS'),
          (R, '(b)  top surface', 'log DOS'),
          (ONLY, '(c)  bottom surface only  '
                 '$\\rho_{\\rm surf}-\\rho_{\\rm bulk}N_{\\rm top}/N_{\\rm dim}$',
           '$\\log_{10}$ surface DOS')]
fig, axs = plt.subplots(1, 3, figsize=(17.8, 6.0), sharey=True)
for ax, (V, tt, cl) in zip(axs, panels):
    lo, hi = np.percentile(V, [2, 99.3])
    im = ax.pcolormesh(kb, om, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    ax.axhline(W0, color='white', lw=1.7, ls='--', alpha=0.9, path_effects=STK)
    for s in (+1, -1):
        ax.plot(s * KBN, W0, 'o', ms=13, mfc='none', mec='white', mew=3.2, zorder=6)
        ax.plot(s * KBN, W0, 'o', ms=13, mfc='none', mec='k', mew=1.3, zorder=7)
        ax.axvline(s * KBN, color='white', lw=1.1, ls=':', alpha=0.75)
    ax.annotate('$\\chi=+1$', xy=(+KBN, W0), xytext=(30, -34), textcoords='offset points',
                color='white', fontsize=13, fontweight='bold', ha='center',
                zorder=8, path_effects=STK)
    ax.annotate('$\\chi=-1$', xy=(-KBN, W0), xytext=(-30, -34), textcoords='offset points',
                color='white', fontsize=13, fontweight='bold', ha='center',
                zorder=8, path_effects=STK)
    ax.set_xlim(KB0, KB1)
    ax.set_xlabel('$k_b$  (reduced)   at  $k_a$ = %.5f' % KA)
    ax.set_title(tt, fontsize=14.5, pad=8)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(cl, fontsize=11.5)
axs[0].set_ylabel('Frequency (THz)')
axs[0].text(0.025, 0.975, '$\\omega_{\\rm Weyl}$ = %.4f THz' % W0, color='white',
            transform=axs[0].transAxes, va='top', fontsize=12.5, path_effects=STK)
fig.suptitle('Surface spectrum on the line through both Weyl projections  '
             '($k_a$ = %.5f,  surface $\\perp$ $c$,  $\\eta$ = 0.005 THz)\n'
             'the two bulk cones cross at $\\omega_{\\rm Weyl}$ exactly at $k_b = \\pm0.07085$ '
             '— the slab Green function reproduces the bulk node position' % KA,
             y=1.035, fontsize=15)
fig.savefig('figs/fig4_slab_arc.png')
print('fig4 저장')
