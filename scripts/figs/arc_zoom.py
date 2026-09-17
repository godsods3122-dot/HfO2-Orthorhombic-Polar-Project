#!/usr/bin/env python3
"""Fine-grid arc map + matching bulk gap mask on the SAME grid.

Answers one question: does the in-gap ridge seen in the m0.8 raw slab terminate
at the projected Weyl nodes, or does it run through them into the continuum?
The coarse 61x61 map (dk=0.0167) cannot resolve that; the node sits between grid
rows.  Here dk=0.005 near the node.

args: SRC HR OUT NORMAL NSLAB W0 ETA  KA0 KA1 NKA  KB0 KB1 NKB  NKC
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/HfO2-Orthorhombic-Polar-Project/scripts/figs')
from slab_asr import read_hr
from phonopy.interface.calculator import read_crystal_structure

SRC, HR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NORMAL, NSLAB = int(sys.argv[4]), int(sys.argv[5])
W0, ETA = float(sys.argv[6]), float(sys.argv[7])
KA0, KA1, NKA = float(sys.argv[8]), float(sys.argv[9]), int(sys.argv[10])
KB0, KB1, NKB = float(sys.argv[11]), float(sys.argv[12]), int(sys.argv[13])
NKC = int(sys.argv[14])

u, _ = read_crystal_structure(filename=SRC + '/POSCAR', interface_mode='vasp')
sqm = np.repeat(np.sqrt(np.array(u.masses)), 3)
R, H, deg = read_hr(HR)
nw = H.shape[1]
inplane = [i for i in range(3) if i != NORMAL]

groups = {}
Hn = H / deg[:, None, None]
for ir in range(len(R)):
    d = int(R[ir, NORMAL])
    groups.setdefault(d, []).append((R[ir, inplane[0]], R[ir, inplane[1]], Hn[ir]))
for d in groups:                              # pack for speed
    r1 = np.array([g[0] for g in groups[d]], float)
    r2 = np.array([g[1] for g in groups[d]], float)
    hh = np.array([g[2] for g in groups[d]])
    groups[d] = (r1, r2, hh)

M = np.tile(sqm, NSLAB); MM = np.outer(M, M)
na = NSLAB * nw // 3
MMb = np.outer(sqm, sqm)

def blocks_at(ka, kb):
    out = {}
    for d, (r1, r2, hh) in groups.items():
        ph = np.exp(2j * np.pi * (ka * r1 + kb * r2))
        out[d] = np.tensordot(ph, hh, axes=(0, 0))
    return out

def slab_D(bl):
    D = np.zeros((NSLAB * nw, NSLAB * nw), complex)
    for d, B in bl.items():
        for n in range(NSLAB):
            m = n + d
            if 0 <= m < NSLAB:
                D[n*nw:(n+1)*nw, m*nw:(m+1)*nw] += B
    return D

def enforce_asr(D):
    P = D * MM
    for i in range(na):
        s = slice(3*i, 3*i+3)
        blk = P[s, :].reshape(3, na, 3)
        P[s, s] = -(blk.sum(axis=1) - blk[:, i, :])
    return P / MM

def bulk_gapped(bl):
    """True if no bulk band crosses W0 for any k_normal on this (ka,kb) rod."""
    lo, hi = [], []
    for t in np.linspace(0, 1, NKC, endpoint=False):
        D = np.zeros((nw, nw), complex)
        for d, B in bl.items():
            D += B * np.exp(2j * np.pi * t * d)
        D = (D + D.conj().T) / 2
        w = np.linalg.eigvalsh(D)
        f = np.sign(w) * np.sqrt(np.abs(w))
        lo.append(f)
    F = np.array(lo)                       # (NKC, nw) band-resolved
    F = np.sort(F, axis=1)
    return not np.any((F.min(axis=0) <= W0) & (F.max(axis=0) >= W0))

kas = np.linspace(KA0, KA1, NKA)
kbs = np.linspace(KB0, KB1, NKB)
res = {k: np.zeros((NKA, NKB)) for k in ('raw', 'asr', 'raw_d', 'asr_d')}
gap = np.zeros((NKA, NKB), bool)
for i, ka in enumerate(kas):
    for j, kb in enumerate(kbs):
        bl = blocks_at(ka, kb)
        gap[i, j] = bulk_gapped(bl)
        D0 = slab_D(bl)
        for tag, D in (('raw', D0), ('asr', enforce_asr(D0))):
            Dh = (D + D.conj().T) / 2
            w, v = np.linalg.eigh(Dh)
            f = np.sign(w) * np.sqrt(np.abs(w))
            top = (np.abs(v[:nw, :])**2).sum(axis=0)
            bot = (np.abs(v[-nw:, :])**2).sum(axis=0)
            L = ETA / ((f - W0)**2 + ETA**2)
            res[tag][i, j] = (top * L).sum()
            res[tag+'_d'][i, j] = ((top - bot) * L).sum()
    print('  ka row %d/%d  (ka=%+.4f)' % (i+1, NKA, ka), flush=True)

np.savez(OUT, kas=kas, kbs=kbs, gap=gap, **res)
print('saved', OUT, ' gap fraction %.3f' % gap.mean())
