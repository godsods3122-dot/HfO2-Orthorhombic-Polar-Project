#!/usr/bin/env python3
"""Arc map restricted to the projected-bulk-gap region.

Outside the gap the surface signal is a resonance mixed with the bulk continuum
and cannot be read as an arc; inside the gap a genuine arc would show as a sharp
line. So mask to the gap and set the colour scale from that region alone.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

npz, maskf, title, out = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
A, B = float(sys.argv[5]), float(sys.argv[6])
d = np.load(npz); ks = d['ks']
gap = np.load(maskf)                      # True = gapped, indexed [ka, kb]

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4))
for ax, tag, lab in zip(axes, ['raw_d', 'asr_d'],
                        ['raw truncation (ASR violated)',
                         'ASR enforced (surface self-term rebuilt)']):
    Z = np.where(gap, d[tag], np.nan).T    # keep gap only
    v = np.nanpercentile(np.abs(Z), 99)
    cmap = plt.get_cmap('RdBu_r').copy(); cmap.set_bad('0.85')
    im = ax.pcolormesh(ks, ks, np.clip(Z, -v, v), cmap=cmap,
                       vmin=-v, vmax=v, shading='auto')
    ax.contour(ks, ks, gap.T.astype(float), levels=[0.5], colors='k',
               linewidths=0.9, alpha=0.6)
    for (a, b, c) in [(A, B, '+'), (-A, -B, '+'), (-A, B, '-'), (A, -B, '-')]:
        ax.plot(a, b, 'o', mfc='none', mec='k', mew=2.0, ms=13)
        ax.text(a, b + 0.05, c, ha='center', va='bottom', fontsize=13, fontweight='bold')
    ax.set_xlabel('$k_a$ (reduced)'); ax.set_ylabel('$k_b$ (reduced)')
    ax.set_title(lab, fontsize=10)
    ax.set_aspect('equal'); ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
    fig.colorbar(im, ax=ax, label='top - bottom surface weight')
fig.suptitle(title + '   (grey = projected bulk continuum, masked out)', fontsize=11)
fig.tight_layout()
fig.savefig(out, dpi=160)
print('wrote', out)
