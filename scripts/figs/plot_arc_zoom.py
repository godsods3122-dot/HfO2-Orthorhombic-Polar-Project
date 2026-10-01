#!/usr/bin/env python3
"""arc_zoom.py 결과(미세 격자, 마스크 동봉) 그림.

arc_zoom.py 는 아크 맵과 간극 마스크를 **같은 격자에서 같은 인덱스 순서로** 뽑으므로
전치 사고(SUMMARY 7절)가 원리적으로 일어나지 않는다.

사용: python3 plot_arc_zoom.py ZOOM.npz TITLE OUT.png N1 N2
      (N1, N2 = 노드 사영 좌표, 배열 인덱스 순서)
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

npz, title, out = sys.argv[1], sys.argv[2], sys.argv[3]
A, B = float(sys.argv[4]), float(sys.argv[5])
d = np.load(npz)
x, y, gap = d['kas'], d['kbs'], d['gap']

fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4))
for ax, tag, lab in zip(axes, ['raw_d', 'asr_d'],
                        ['raw truncation (ASR violated)',
                         'ASR enforced (surface self-term rebuilt)']):
    Z = np.where(gap, d[tag], np.nan).T
    v = np.nanpercentile(np.abs(Z), 99)
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('0.85')
    im = ax.pcolormesh(x, y, np.clip(Z, -v, v), cmap=cmap, vmin=-v, vmax=v, shading='auto')
    ax.contour(x, y, gap.T.astype(float), levels=[0.5], colors='k', linewidths=0.9, alpha=0.6)
    for (a, b, c) in [(A, B, '+'), (A, -B, '-')]:
        ax.plot(a, b, 'o', mfc='none', mec='k', mew=2.0, ms=13)
        ax.text(a, b + 0.02, c, ha='center', va='bottom', fontsize=13, fontweight='bold')
    ax.set_xlabel('$k_1$ (reduced, array index 0)')
    ax.set_ylabel('$k_2$ (reduced, array index 1)')
    ax.set_title(lab, fontsize=10)
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax, label='top - bottom surface weight')
fig.suptitle(title + '   (grey = projected bulk continuum)', fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=160)
print('wrote', out)
