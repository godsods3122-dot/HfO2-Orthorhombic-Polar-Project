#!/usr/bin/env python3
"""아크 데이터와 간극 마스크의 축 규약이 서로 맞는지 검사한다.

이걸 안 해서 SUMMARY 7절이 두 번 틀렸다. m0.8 의 마스크가 아크 배열에 대해 전치돼
있었고, 그 결과 있지도 않은 "초승달"이 그림에 나타났다.

두 가지를 독립적으로 확인한다.
  1) 아크 배열의 법선축과 인덱스 순서 — 몇 개 격자점에서 raw_d 를 재계산해 대조.
     맞는 조합에서만 차이가 정확히 0 이 된다.
  2) 마스크의 인덱스 순서 — 같은 법선으로 마스크를 재계산해 그대로/전치 일치율 비교.
     마스크가 대칭이면 이 검사가 무의미하므로 저장본의 자기전치 일치율도 함께 낸다.

사용: python3 arc_axis_check.py SRC HR ARC.npz [GAPMASK.npy] W0 [NSLAB] [ETA]
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from slab_asr import read_hr
from phonopy.interface.calculator import read_crystal_structure

SRC, HRF, ARC = sys.argv[1], sys.argv[2], sys.argv[3]
MASK = sys.argv[4] if not sys.argv[4].replace('.', '').isdigit() else None
rest = sys.argv[5:] if MASK else sys.argv[4:]
W0 = float(rest[0])
NSLAB = int(rest[1]) if len(rest) > 1 else 16
ETA = float(rest[2]) if len(rest) > 2 else 0.02

u, _ = read_crystal_structure(filename=SRC + '/POSCAR', interface_mode='vasp')
R, H, deg = read_hr(HRF)
nw = H.shape[1]
Hn = H / deg[:, None, None]
arc = np.load(ARC)
ks = arc['ks']
NG = len(ks)


def group(normal):
    ip = [i for i in range(3) if i != normal]
    g = {}
    for ir in range(len(R)):
        d = int(R[ir, normal])
        g.setdefault(d, [[], [], []])
        g[d][0].append(R[ir, ip[0]]); g[d][1].append(R[ir, ip[1]]); g[d][2].append(Hn[ir])
    return {d: (np.array(v[0], float), np.array(v[1], float), np.array(v[2]))
            for d, v in g.items()}


def raw_d_at(G, x, y):
    D = np.zeros((NSLAB * nw, NSLAB * nw), complex)
    for d, (r1, r2, hh) in G.items():
        B = np.tensordot(np.exp(2j * np.pi * (x * r1 + y * r2)), hh, axes=(0, 0))
        for n in range(NSLAB):
            m = n + d
            if 0 <= m < NSLAB:
                D[n * nw:(n + 1) * nw, m * nw:(m + 1) * nw] += B
    D = (D + D.conj().T) / 2
    w, v = np.linalg.eigh(D)
    f = np.sign(w) * np.sqrt(np.abs(w))
    top = (np.abs(v[:nw, :]) ** 2).sum(0)
    bot = (np.abs(v[-nw:, :]) ** 2).sum(0)
    L = ETA / ((f - W0) ** 2 + ETA ** 2)
    return ((top - bot) * L).sum()


def gapmask(G, nkc=41):
    gap = np.zeros((NG, NG), bool)
    ts = np.linspace(0, 1, nkc, endpoint=False)
    for i, ka in enumerate(ks):
        for j, kb in enumerate(ks):
            F = []
            for t in ts:
                D = np.zeros((nw, nw), complex)
                for d, (r1, r2, hh) in G.items():
                    D += np.tensordot(np.exp(2j * np.pi * (ka * r1 + kb * r2 + t * d)),
                                      hh, axes=(0, 0))
                D = (D + D.conj().T) / 2
                w = np.linalg.eigvalsh(D)
                F.append(np.sort(np.sign(w) * np.sqrt(np.abs(w))))
            F = np.array(F)
            gap[i, j] = not np.any((F.min(0) <= W0) & (F.max(0) >= W0))
    return gap


tests = [(NG // 3, NG // 2), (NG // 4, 2 * NG // 3), (NG // 2, NG // 2),
         (2 * NG // 3, NG // 4), (3 * NG // 5, 5 * NG // 11)]
print('1) 아크 배열의 법선축·인덱스 순서   (차이 0 인 조합이 정답)')
best = None
for normal in (0, 1, 2):
    G = group(normal)
    r = [raw_d_at(G, ks[i], ks[j]) for i, j in tests]
    e_ij = max(abs(r[t] - arc['raw_d'][i, j]) for t, (i, j) in enumerate(tests))
    e_ji = max(abs(r[t] - arc['raw_d'][j, i]) for t, (i, j) in enumerate(tests))
    print('   NORMAL=%d   [i,j] 차이 %9.4f    [j,i] 차이 %9.4f' % (normal, e_ij, e_ji))
    if best is None or min(e_ij, e_ji) < best[1]:
        best = (normal, min(e_ij, e_ji), e_ij <= e_ji)
print('   -> NORMAL=%d, 인덱스 %s' % (best[0], '[i,j] 그대로' if best[2] else '[j,i] 전치'))

if MASK:
    ref = np.load(MASK)
    g = gapmask(group(best[0]))
    print()
    print('2) 마스크의 인덱스 순서')
    print('   그대로 일치 %.4f   전치하면 일치 %.4f   (저장본 자기전치 일치 %.4f)'
          % ((g == ref).mean(), (g.T == ref).mean(), (ref == ref.T).mean()))
    if (ref == ref.T).mean() > 0.98:
        print('   -> 마스크가 대칭이라 판정 불가')
    elif (g.T == ref).mean() > (g == ref).mean():
        print('   -> 마스크가 **전치돼 있다**. 그대로 쓰면 안 된다.')
    else:
        print('   -> 정상')
