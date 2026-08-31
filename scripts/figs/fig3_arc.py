#!/usr/bin/env python3
"""Fig 3: ω = 10.0869 THz 등주파수 표면 스펙트럼 (bottom / top / bulk + 아크 분리).

Simphony SlabArc_calc 의 arc.dat_l / arc.dat_r / arc.dat_bulk 를 읽는다.
표면 법선은 c(편극축)이라 표면 BZ 좌표가 (k_a, k_b) 그대로다.

이 물질은 노드 주파수에서 사영 벌크 연속체가 표면 BZ 전체를 덮는다
(phonopy 로 확인: 창 안 100 %).  그래서 페르미 아크는 간극 속 고립 상태가
아니라 **공명**이고, 벌크 배경이 큰 raw DOS 만 보면 안 보인다.
위/아래 표면 차 `dos_l − dos_r` 가 아크만 남긴다 (벌크는 상쇄된다).

겹쳐 그린 곡선은 phonopy 로 따로 구한 **band 17 / band 18 의 사영 연속체
경계**다.  두 경계가 교차하는 자리가 정확히 바일 사영점이고, 아크는 거기서
끝난다.  이게 아크임을 보이는 결정적 증거다.

⚠️ Simphony `fermiarc.f90` 은 KPLANE_SLAB 의 두 벡터를 정수로 올림한다
(`ceiling`/`floor`, QPI 용). 그래서 좁은 창을 줘도 실제로는 K2D_start 에서
시작하는 full BZ 한 칸을 계산한다.  파일의 kx, ky 폭이 곧 |b1|, |b2| 이므로
그걸로 나누면 reduced 좌표가 나온다.  여기서는 [-0.5,0.5) 로 되감는다.
(`patch/apply_fermiarc_kplane_fix.py` 참조)
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

D = sys.argv[1] if len(sys.argv) > 1 else 'work/slabarc'
XL, YL = (-0.33, 0.33), (-0.23, 0.23)
NODES = [(+0.1464927, +0.0708493, +1), (+0.1464927, -0.0708493, -1),
         (-0.1464927, +0.0708493, -1), (-0.1464927, -0.0708493, +1)]
EARC = 10.0869
PROJ = 'figs/proj_continuum_1718.npz'


def load(name):
    d = np.loadtxt(os.path.join(D, name))
    kx, ky, v = d[:, 0], d[:, 1], d[:, 2]
    n_in = int(np.argmax(np.diff(ky) < 0) + 1)          # ky 가 안쪽 루프
    n_out = len(v) // n_in
    V = v.reshape(n_out, n_in)
    b1, b2 = kx.max() - kx.min(), ky.max() - ky.min()   # = |b1|, |b2| (full BZ)
    a = kx.reshape(n_out, n_in)[:, 0] / b1
    b = ky.reshape(n_out, n_in)[0] / b2
    a, b, V = a[:-1], b[:-1], V[:-1, :-1]               # 마지막은 주기 복제
    a, b = (a + 0.5) % 1.0 - 0.5, (b + 0.5) % 1.0 - 0.5
    ia, ib = np.argsort(a), np.argsort(b)
    return a[ia], b[ib], V[np.ix_(ia, ib)]


def projected():
    """band 17 / 18 의 사영 연속체 (k_c 로 사영). 없으면 phonopy 로 만든다."""
    if os.path.exists(PROJ):
        z = np.load(PROJ)
        return z['A'], z['B'], z['in17'], z['in18']
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))))
    from weyl_scan import get_ph
    ph = get_ph('source/parent_pristine')
    NA, NB, NC = 137, 97, 41
    A = np.linspace(-0.34, 0.34, NA); B = np.linspace(-0.24, 0.24, NB)
    C = np.linspace(-0.5, 0.5, NC)
    ph.run_qpoints([[x, y, z] for x in A for y in B for z in C])
    f = ph.qpoints.frequencies.reshape(NA, NB, NC, -1)
    w0 = 10.086936
    i17 = (f[..., 16].min(2) <= w0) & (f[..., 16].max(2) >= w0)
    i18 = (f[..., 17].min(2) <= w0) & (f[..., 17].max(2) >= w0)
    np.savez(PROJ, A=A, B=B, in17=i17, in18=i18)
    return A, B, i17, i18


a, b, L = load('arc.dat_l')
_, _, R = load('arc.dat_r')
_, _, BK = load('arc.dat_bulk')
PA, PB, in17, in18 = projected()

panels = [(L, '(a)  bottom surface', 'log DOS'),
          (R, '(b)  top surface', 'log DOS'),
          (BK, '(c)  bulk projection', 'log DOS'),
          (L - R, '(d)  bottom $-$ top  (Fermi arc isolated)', '$\\log(\\rho_{\\rm bot}/\\rho_{\\rm top})$')]
fig, axs = plt.subplots(1, 4, figsize=(23.0, 5.9), sharex=True, sharey=True)
for ax, (V, tt, cl) in zip(axs, panels):
    m = (a >= XL[0]) & (a <= XL[1]); n = (b >= YL[0]) & (b <= YL[1])
    lo, hi = np.percentile(V[np.ix_(m, n)], [2, 99])
    im = ax.pcolormesh(a, b, V.T, cmap='jet', vmin=lo, vmax=hi,
                       shading='gouraud', rasterized=True)
    c17 = ax.contour(PA, PB, in17.T.astype(float), levels=[0.5],
                     colors='white', linewidths=1.9)
    c18 = ax.contour(PA, PB, in18.T.astype(float), levels=[0.5],
                     colors='black', linewidths=1.9, linestyles='dashed')
    for cc, fg in ((c17, 'black'), (c18, 'white')):
        for col in (cc.collections if hasattr(cc, 'collections') else [cc]):
            col.set_path_effects([pe.withStroke(linewidth=3.6, foreground=fg)])
    for ka, kb, c in NODES:
        ax.plot(ka, kb, 'o', ms=13, mfc='none', mec='white', mew=3.4, zorder=6)
        ax.plot(ka, kb, 'o', ms=13, mfc='none', mec='k', mew=1.4, zorder=7)
        ax.annotate('$\\chi=%+d$' % c, xy=(ka, kb),
                    xytext=(34 if ka > 0 else -34, 20 if kb > 0 else -20),
                    textcoords='offset points', color='white',
                    ha='center', va='center', fontsize=13.5, fontweight='bold', zorder=8,
                    path_effects=[pe.withStroke(linewidth=3.2, foreground='black')])
    ax.set_xlabel('$k_a$  (reduced)')
    ax.set_title(tt, fontsize=14.5, pad=8)
    ax.set_xlim(*XL); ax.set_ylim(*YL)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(cl, fontsize=11.5)
axs[0].set_ylabel('$k_b$  (reduced)')

from matplotlib.lines import Line2D
STK = [pe.withStroke(linewidth=3.6, foreground='black')]
axs[0].legend(handles=[Line2D([], [], color='white', lw=2.0, path_effects=STK,
                              label='band 17 projected edge'),
                       Line2D([], [], color='black', lw=2.0, ls=(0, (5, 3)),
                              path_effects=[pe.withStroke(linewidth=3.6, foreground='white')],
                              label='band 18 projected edge')],
              loc='upper center', bbox_to_anchor=(0.55, -0.155), ncol=2,
              frameon=True, facecolor='#cfcfcf', edgecolor='#888888',
              framealpha=1.0, fontsize=12)
axs[3].text(0.5, 0.955, 'Fermi arc', transform=axs[3].transAxes, ha='center', va='top',
            color='white', fontsize=13.5, fontweight='bold', zorder=9, path_effects=STK)
for sgn in (+1, -1):
    axs[3].annotate('', xy=(sgn * 0.243, 0.163), xytext=(sgn * 0.055, 0.198),
                    arrowprops=dict(arrowstyle='-|>', color='white', lw=2.0,
                                    shrinkB=3, path_effects=STK), zorder=9)

fig.suptitle('Iso-frequency surface spectrum at the Weyl energy  $\\omega$ = %.4f THz'
             '   (surface $\\perp$ polar axis $c$,  $N_p$ = 2,  $\\eta$ = 0.008 THz)\n'
             'the two projected band edges cross exactly at the four Weyl projections — '
             'the arc terminates there' % EARC, y=1.06, fontsize=15)
fig.savefig('figs/fig3_fermi_arc.png')
print('fig3 저장')
