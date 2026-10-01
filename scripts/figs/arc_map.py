#!/usr/bin/env python3
"""표면 BZ 위의 표면-사영 스펙트럼 세기 (아크 맵).

슬랩을 hr.dat(총 힘상수)에서 직접 잘라 만들고, 두 가지로 계산한다:
  raw : 자른 그대로 (ASR 위반 -> Simphony 현재 동작과 같은 결함)
  asr : 표면 자기항을 sum_j Phi_ij = 0 이 되도록 재구성 (ASR 정확)
둘을 비교해 아크 형태가 설정에 얼마나 좌우되는지 본다.
"""
import sys, numpy as np
sys.path.insert(0, '/tmp/claude-0/-home-user-HfO2-Orthorhombic-Polar-Project/c2573bf3-53de-56be-b0b9-1d949d54780e/scratchpad')
from slab_asr import read_hr
from phonopy.interface.calculator import read_crystal_structure

SRC, HR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NORMAL = int(sys.argv[4]); NSLAB = int(sys.argv[5]); NG = int(sys.argv[6])
W0 = float(sys.argv[7]); ETA = float(sys.argv[8])

u, _ = read_crystal_structure(filename=SRC + '/POSCAR', interface_mode='vasp')
sqm = np.repeat(np.sqrt(np.array(u.masses)), 3)
R, H, deg = read_hr(HR)
nw = H.shape[1]
inplane = [i for i in range(3) if i != NORMAL]

# 법선 R 성분별로, 면내 R 은 유지
groups = {}
for ir in range(len(R)):
    d = int(R[ir, NORMAL])
    groups.setdefault(d, []).append((R[ir, inplane[0]], R[ir, inplane[1]], H[ir] / deg[ir]))

M = np.tile(sqm, NSLAB)
MM = np.outer(M, M)
na = NSLAB * nw // 3


def slab_D(ka, kb):
    D = np.zeros((NSLAB * nw, NSLAB * nw), complex)
    for d, lst in groups.items():
        B = np.zeros((nw, nw), complex)
        for r1, r2, h in lst:
            B += h * np.exp(2j * np.pi * (ka * r1 + kb * r2))
        for n in range(NSLAB):
            m = n + d
            if 0 <= m < NSLAB:
                D[n * nw:(n + 1) * nw, m * nw:(m + 1) * nw] += B
    return D


def enforce_asr(D):
    P = D * MM
    for i in range(na):
        s = slice(3 * i, 3 * i + 3)
        blk = P[s, :].reshape(3, na, 3)
        P[s, s] = -(blk.sum(axis=1) - blk[:, i, :])
    return P / MM


ks = np.linspace(-0.5, 0.5, NG)
surf = max(1, nw // 3)          # 최상단 1 유닛셀 층의 궤도 수
res = {k: np.zeros((NG, NG)) for k in ('raw', 'asr', 'raw_d', 'asr_d')}
for i, ka in enumerate(ks):
    for j, kb in enumerate(ks):
        D0 = slab_D(ka, kb)
        for tag, D in (('raw', D0), ('asr', enforce_asr(D0))):
            Dh = (D + D.conj().T) / 2
            w, v = np.linalg.eigh(Dh)
            f = np.sign(w) * np.sqrt(np.abs(w))
            top = (np.abs(v[:nw, :]) ** 2).sum(axis=0)
            bot = (np.abs(v[-nw:, :]) ** 2).sum(axis=0)
            L = ETA / ((f - W0) ** 2 + ETA ** 2)
            res[tag][i, j] = (top * L).sum()
            res[tag + '_d'][i, j] = ((top - bot) * L).sum()
    if i % 10 == 0:
        print('  k row %d/%d' % (i, NG), flush=True)

np.savez(OUT, ks=ks, **res)
print('saved', OUT)
