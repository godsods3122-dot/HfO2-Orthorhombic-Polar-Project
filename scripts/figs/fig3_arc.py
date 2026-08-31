#!/usr/bin/env python3
"""Fig 3: ω = 10.0869 THz 등주파수 표면 스펙트럼 (bottom / top / bulk / 표면 전용).

Simphony `SlabArc_calc` 의 arc.dat_l / arc.dat_r / arc.dat_bulk 를 읽는다.
표면 법선은 c(편극축)이라 표면 BZ 좌표가 (k_a, k_b) 그대로다.

## 이 그림을 제대로 뽑기 위해 고쳐야 했던 것 세 가지

1. **LO-TO 누락** (`patch/apply_fermiarc_loto_fix.py`).
   `fermiarc.f90` 은 `ham_qlayer2qlayer` 를 조건 없이 불러서 SlabArc 가
   LO-TO 없이 계산됐다.  `surfstat.f90` 은 제대로 갈래를 탄다.
2. **KPLANE_SLAB 무시** (`patch/apply_fermiarc_kplane_fix.py`).
   두 벡터를 정수로 올려(ceiling) 항상 full BZ 를 계산했다.
3. **표면 전용 DOS 정규화** (`patch/apply_surfdos_only_norm_fix.py`).
   `dos_l` 은 표면 단위포 궤도(NtopOrbitals)만, `dos_bulk` 은 주층 전체
   (Ndim = Np*Num_wann)를 합한다.  그냥 빼면 벌크를 Np 배로 빼는 셈이라
   결과가 항상 음수다.  `dos_bulk * NtopOrbitals/Ndim` 로 맞춰야 한다.
   fermiarc 에는 `_only` 출력 자체가 없어 여기서 후처리로 계산한다.

노드 주파수에는 벌크 간극이 없다 (사영 연속체 커버리지 1.000).  그래서
raw DOS (a,b,c) 로는 아크가 벌크에 묻힌다.  (d) 의 표면 전용 값에서
표면 상태 영역이 날카롭게 잘리고, 그 경계가 네 바일 사영점을 지난다.

겹쳐 그린 곡선은 phonopy 로 따로 구한 band 17 / band 18 의 사영 연속체
경계다.  두 경계가 교차하는 자리가 정확히 바일 사영점이다.
"""
import sys, os, glob, re
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import setup
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

DIRS = sys.argv[1:] or sorted(glob.glob('work/arc_loto_*'))
XL, YL = (-0.30, 0.30), (-0.20, 0.20)
NODES = [(+0.1464927, +0.0708493, +1), (+0.1464927, -0.0708493, -1),
         (-0.1464927, +0.0708493, -1), (-0.1464927, -0.0708493, +1)]
EARC = 10.0869
PROJ = 'figs/proj_continuum_1718.npz'
STK = [pe.withStroke(linewidth=3.4, foreground='black')]


def orbital_ratio(d):
    """NtopOrbitals / Ndim — 표면 합과 벌크 합의 궤도 개수 비."""
    s = open(os.path.join(d, 'PN.out')).read()
    nt = int(re.search(r'NtopOrbitals\s+(\d+)', s).group(1))
    nd = int(re.search(r'(?i)\bndim:?\s+(\d+)', s).group(1))
    return nt / nd


def load_one(d, name):
    a = np.loadtxt(os.path.join(d, name))
    kx, ky, v = a[:, 0], a[:, 1], a[:, 2]
    n_in = int(np.argmax(np.diff(ky) < 0) + 1)          # ky 가 안쪽 루프
    n_out = len(v) // n_in
    # KPLANE 패치 후에는 파일 좌표가 지정한 창 그대로다.  reduced 로 되돌리려면
    # PN.out 의 K2D_start / 벡터를 쓴다.
    s = open(os.path.join(d, 'PN.out')).read()
    st = [float(x) for x in re.search(r'K2D_start:\s*(\S+)\s+(\S+)', s).groups()]
    v1 = [float(x) for x in re.search(r'The first vector:\s*(\S+)\s+(\S+)', s).groups()]
    v2 = [float(x) for x in re.search(r'The second vector:\s*(\S+)\s+(\S+)', s).groups()]
    ka = st[0] + v1[0] * np.linspace(0, 1, n_out)
    kb = st[1] + v2[1] * np.linspace(0, 1, n_in)
    return ka, kb, v.reshape(n_out, n_in)


def stitch(name):
    """k_a 로 쪼개 돌린 조각들을 하나로 붙인다 (경계 중복 열은 버린다)."""
    parts = [load_one(d, name) for d in DIRS]
    parts.sort(key=lambda p: p[0][0])
    kb = parts[0][1]
    ka = parts[0][0]; V = parts[0][2]
    for k2, _, V2 in parts[1:]:
        if abs(k2[0] - ka[-1]) < 1e-9:
            k2, V2 = k2[1:], V2[1:]
        ka = np.concatenate([ka, k2]); V = np.concatenate([V, V2], axis=0)
    return ka, kb, V


a, b, L = stitch('arc.dat_l')
_, _, R = stitch('arc.dat_r')
_, _, BK = stitch('arc.dat_bulk')
ratio = orbital_ratio(DIRS[0])
only_l = np.maximum(np.exp(L) - np.exp(BK) * ratio, 0.0)
z = np.load(PROJ)
PA, PB, in17, in18 = z['A'], z['B'], z['in17'], z['in18']

panels = [(L, '(a)  bottom surface', 'log DOS', 'jet'),
          (R, '(b)  top surface', 'log DOS', 'jet'),
          (BK, '(c)  bulk projection', 'log DOS', 'jet'),
          (np.log10(only_l + 1e-3),
           '(d)  bottom surface only  $\\rho_{\\rm surf}-\\rho_{\\rm bulk}\\,N_{\\rm top}/N_{\\rm dim}$',
           '$\\log_{10}$ surface DOS', 'jet')]
fig, axs = plt.subplots(1, 4, figsize=(23.0, 6.2), sharex=True, sharey=True)
for ax, (V, tt, cl, cm) in zip(axs, panels):
    m = (a >= XL[0]) & (a <= XL[1]); n = (b >= YL[0]) & (b <= YL[1])
    lo, hi = np.percentile(V[np.ix_(m, n)], [3, 99.3])
    im = ax.pcolormesh(a, b, V.T, cmap=cm, vmin=lo, vmax=hi,
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
                    xytext=(36 if ka > 0 else -36, 22 if kb > 0 else -22),
                    textcoords='offset points', color='white',
                    ha='center', va='center', fontsize=13, fontweight='bold',
                    zorder=8, path_effects=STK)
    ax.set_xlabel('$k_a$  (reduced)')
    ax.set_title(tt, fontsize=13.5, pad=8)
    ax.set_xlim(*XL); ax.set_ylim(*YL)
    cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
    cb.set_label(cl, fontsize=11.5)
axs[0].set_ylabel('$k_b$  (reduced)')
axs[3].text(0.5, 0.955, 'Fermi arc', transform=axs[3].transAxes, ha='center',
            va='top', color='white', fontsize=13.5, fontweight='bold',
            zorder=9, path_effects=STK)

fig.legend(handles=[Line2D([], [], color='white', lw=2.0, path_effects=STK,
                           label='band 17 projected edge'),
                    Line2D([], [], color='black', lw=2.0, ls='dashed',
                           path_effects=[pe.withStroke(linewidth=3.6, foreground='white')],
                           label='band 18 projected edge')],
           loc='lower center', bbox_to_anchor=(0.5, -0.035), ncol=2,
           frameon=True, facecolor='#cfcfcf', edgecolor='#888888',
           framealpha=1.0, fontsize=12.5)
fig.suptitle('Iso-frequency surface spectrum at the Weyl energy  $\\omega$ = %.4f THz'
             '   (surface $\\perp$ polar axis $c$,  LO-TO on,  $N_p$ = 2,  $\\eta$ = 0.004 THz)\n'
             'the two projected band edges cross exactly at the four Weyl projections '
             '— the surface-state region is bounded by them and terminates there' % EARC,
             y=1.055, fontsize=15)
fig.savefig('figs/fig3_fermi_arc.png')
print('fig3 저장.  격자 %d x %d,  간격 %.5f x %.5f,  N_top/N_dim = %.3f'
      % (len(a), len(b), a[1] - a[0], b[1] - b[0], ratio))
