#!/usr/bin/env python3
"""발표용 — −0.8 % / 무변형 / +1 % 밴드 중첩. 미러 배치 세 점만.

compressive 빨강, tensile 파랑, 무변형 회색. 전부 반투명.
x 축은 세 구조가 같은 환산 경로를 쓰므로 무변형 구조의 경로 길이로 통일한다.
안쪽 Γ 는 LO-TO 방향 의존 불연속이라 구조마다 NaN 으로 끊는다.

출력: figs/fig_overlay_ppt.png  (데이터 캐시: figs/overlay_bands.npz)
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from axes import roles, to_native, HS, PATH
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

NK = 120
CASES = [('$-0.8$ %  compressive', 'm1_mirror',       '#c0392b'),
         ('unstrained',            'pristine_mirror', '#4a4a4a'),
         ('$+1$ %  tensile',       'p1_mirror',       '#1f5fd0')]
CACHE = 'figs/overlay_bands.npz'

if not os.path.exists(CACHE):
    from weyl_scan import get_ph
    out = {}
    for _, src, _ in CASES:
        D = 'source/' + src
        r2n, _, _ = roles(D); ph = get_ph(D)
        rec = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi
        Q, xs, x0, seg = [], [], 0.0, []
        for a, b in zip(PATH[:-1], PATH[1:]):
            ka, kb = np.array(HS[a], float), np.array(HS[b], float)
            L = np.linalg.norm((np.array(to_native(kb, r2n)) - np.array(to_native(ka, r2n))) @ rec)
            t = np.linspace(0, 1, NK, endpoint=False)
            Q += [to_native(ka + tt * (kb - ka), r2n) for tt in t]
            xs.append(x0 + t * L); seg.append(x0); x0 += L
        Q.append(to_native(np.array(HS[PATH[-1]], float), r2n))
        F = []
        for i in range(0, len(Q), 20000):
            ph.run_qpoints([list(q) for q in Q[i:i+20000]])
            F.append(ph.qpoints.frequencies)
        out['E_' + src] = np.concatenate(F).T
        out['x_' + src] = np.append(np.concatenate(xs), x0)
        out['t_' + src] = np.array(seg + [x0])
        print('  %s 계산 완료' % src, flush=True)
    np.savez(CACHE, **out)
d = np.load(CACHE)

REF = 'pristine_mirror'
x, ticks = d['x_' + REF], d['t_' + REF]
LAB = ['Γ' if p == 'G' else p for p in PATH]
GB = 4 * NK

plt.rcParams.update({
    'font.size': 28, 'axes.labelsize': 48, 'legend.fontsize': 27,
    'ytick.labelsize': 34, 'xtick.labelsize': 60,
    'axes.linewidth': 2.4, 'xtick.major.width': 2.4, 'ytick.major.width': 2.4,
    'xtick.major.size': 9, 'ytick.major.size': 9,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'ytick.right': True,
    'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

xp = np.insert(x, GB, np.nan)
fig, ax = plt.subplots(figsize=(15.0, 9.2))
for lab, src, col in CASES:
    E = np.insert(d['E_' + src], GB, np.nan, axis=1)
    for b in range(E.shape[0]):
        if b not in (16, 17):
            ax.plot(xp, E[b], color=col, lw=1.6, alpha=0.22, zorder=2)
    ax.plot(xp, E[16], color=col, lw=4.4, alpha=0.75, zorder=4)
    ax.plot(xp, E[17], color=col, lw=4.4, alpha=0.75, zorder=4)
for t in ticks[1:-1]:
    ax.axvline(t, color='#999999', lw=1.4, alpha=0.55, zorder=1)

ax.set_xticks(ticks); ax.set_xticklabels(LAB, fontsize=60)
ax.set_xlim(x[0], x[-1]); ax.set_ylim(8.8, 12.6)
ax.set_ylabel('Frequency (THz)', labelpad=12)
ax.tick_params(axis='x', pad=10)
ax.legend(handles=[Line2D([], [], color=c, lw=7, alpha=0.8, label=l) for l, _, c in CASES],
          loc='upper center', bbox_to_anchor=(0.5, -0.135), ncol=3,
          frameon=False, handlelength=1.5, columnspacing=1.6)
fig.savefig('figs/fig_overlay_ppt.png')
print('wrote figs/fig_overlay_ppt.png')
