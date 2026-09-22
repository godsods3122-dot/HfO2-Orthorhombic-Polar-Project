#!/usr/bin/env python3
"""발표용 — Γ-X-S 구간만 떼어낸 밴드 중첩, 고해상도.

fig_overlay_ppt.py 와 같은 내용이지만 구간당 400점으로 다시 계산한다
(전 경로판은 120점). 이 구간에는 안쪽 Γ 가 없으므로 끊을 필요도 없다.

압축(−) 파랑 / 무변형 회색 / 인장(+) 빨강, 전부 반투명.
출력: figs/fig_overlay_gxs_ppt.png  (캐시 figs/overlay_gxs.npz)

색 규약: **양수 빨강 / 음수 파랑**. 여기서는 strain 부호에 적용한다 —
인장(+) 빨강, 압축(−) 파랑. 추이 그림(fig8_trend_ppt.py)의 chirality 색도 같은 규칙이다.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from axes import roles, to_native, HS
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

NK = 400
SEGS = [('G', 'X'), ('X', 'S')]
LAB = ['Γ', 'X', 'S']
CASES = [('$-0.8$ %  compressive', 'm1_mirror',       '#1f5fd0'),
         ('unstrained',            'pristine_mirror', '#4a4a4a'),
         ('$+1$ %  tensile',       'p1_mirror',       '#c0392b')]
CACHE = 'figs/overlay_gxs.npz'

if not os.path.exists(CACHE):
    from weyl_scan import get_ph
    out = {}
    for _, src, _ in CASES:
        D = 'source/' + src
        r2n, _, _ = roles(D); ph = get_ph(D)
        rec = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi
        Q, xs, x0, seg = [], [], 0.0, []
        for a, b in SEGS:
            ka, kb = np.array(HS[a], float), np.array(HS[b], float)
            L = np.linalg.norm((np.array(to_native(kb, r2n)) - np.array(to_native(ka, r2n))) @ rec)
            t = np.linspace(0, 1, NK, endpoint=False)
            Q += [to_native(ka + tt * (kb - ka), r2n) for tt in t]
            xs.append(x0 + t * L); seg.append(x0); x0 += L
        Q.append(to_native(np.array(HS[SEGS[-1][1]], float), r2n))
        F = []
        for i in range(0, len(Q), 20000):
            ph.run_qpoints([list(q) for q in Q[i:i+20000]])
            F.append(ph.qpoints.frequencies)
        out['E_' + src] = np.concatenate(F).T
        out['x_' + src] = np.append(np.concatenate(xs), x0)
        out['t_' + src] = np.array(seg + [x0])
        print('  %s (%d점) 완료' % (src, len(Q)), flush=True)
    np.savez(CACHE, **out)
d = np.load(CACHE)

REF = 'pristine_mirror'
x, ticks = d['x_' + REF], d['t_' + REF]

plt.rcParams.update({
    'font.size': 28, 'axes.labelsize': 48, 'legend.fontsize': 27,
    'ytick.labelsize': 34, 'xtick.labelsize': 60,
    'axes.linewidth': 2.4, 'xtick.major.width': 2.4, 'ytick.major.width': 2.4,
    'xtick.major.size': 9, 'ytick.major.size': 9,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'ytick.right': True,
    'savefig.dpi': 400, 'savefig.bbox': 'tight',
    'font.family': 'DejaVu Sans', 'mathtext.fontset': 'dejavusans',
})

fig, ax = plt.subplots(figsize=(13.0, 9.6))
for lab, src, col in CASES:
    E = d['E_' + src]
    for b in range(E.shape[0]):
        if b not in (16, 17):
            ax.plot(x, E[b], color=col, lw=1.8, alpha=0.22, zorder=2)
    ax.plot(x, E[16], color=col, lw=5.0, alpha=0.75, zorder=4)
    ax.plot(x, E[17], color=col, lw=5.0, alpha=0.75, zorder=4)
for t in ticks[1:-1]:
    ax.axvline(t, color='#999999', lw=1.6, alpha=0.55, zorder=1)

ax.set_xticks(ticks); ax.set_xticklabels(LAB, fontsize=60)
ax.set_xlim(x[0], x[-1]); ax.set_ylim(9.2, 12.4)
ax.set_ylabel('Frequency (THz)', labelpad=12)
ax.tick_params(axis='x', pad=10)
ax.legend(handles=[Line2D([], [], color=c, lw=7, alpha=0.8, label=l) for l, _, c in CASES],
          loc='upper center', bbox_to_anchor=(0.5, -0.125), ncol=3,
          frameon=False, handlelength=1.5, columnspacing=1.6)
fig.savefig('figs/fig_overlay_gxs_ppt.png')
print('wrote figs/fig_overlay_gxs_ppt.png   (구간당 %d점)' % NK)
