#!/usr/bin/env python3
"""Surface-selective arc map: (top - bottom) layer weight, raw vs ASR-enforced.

The plain top-layer weight is dominated by the projected bulk continuum; taking the
difference between the two surfaces cancels bulk-like states and leaves the
surface-localised ones.
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

npz, title, out = sys.argv[1], sys.argv[2], sys.argv[3]
A, B = float(sys.argv[4]), float(sys.argv[5])   # node projections (+-A, +-B)
d = np.load(npz)
ks = d['ks']

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.4))
for ax, tag, lab in zip(axes, ['raw_d', 'asr_d'],
                        ['raw truncation (ASR violated)',
                         'ASR enforced (surface self-term rebuilt)']):
    Z = d[tag].T
    v = np.percentile(np.abs(Z), 99.5)
    im = ax.pcolormesh(ks, ks, np.clip(Z, -v, v), cmap='RdBu_r',
                       vmin=-v, vmax=v, shading='auto')
    # node projections: chirality +1 at (+,+),(-,-) and -1 at (-,+),(+,-)
    for (a, b, c) in [(A, B, '+'), (-A, -B, '+'), (-A, B, '-'), (A, -B, '-')]:
        ax.plot(a, b, 'o', mfc='none', mec='k', mew=2.0, ms=13)
        ax.text(a, b + 0.045, c, ha='center', va='bottom', fontsize=12, fontweight='bold')
    # the two candidate arc paths (both shown to be bulk-covered)
    ax.plot([-A, A], [B, B], 'k--', lw=1.2, alpha=0.7)
    ax.plot([A, A], [-B, B], 'k--', lw=1.2, alpha=0.7)
    ax.set_xlabel('$k_a$ (reduced)')
    ax.set_ylabel('$k_b$ (reduced)')
    ax.set_title(lab, fontsize=10)
    ax.set_aspect('equal')
    ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
    fig.colorbar(im, ax=ax, label='top - bottom surface weight')
fig.suptitle(title, fontsize=12)
fig.tight_layout()
fig.savefig(out, dpi=160)
print('wrote', out)
