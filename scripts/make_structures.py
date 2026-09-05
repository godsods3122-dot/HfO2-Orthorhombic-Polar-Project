#!/usr/bin/env python3
"""
make_structures.py — 두 갈래 검증에 필요한 POSCAR 를 생성한다 (POSCAR 만, 계산 없음).

(1) 극성 진폭 손잡이  r0 + lambda*u
    lambda = 1-2x (x = mirror 결함 분율). Pca2_1 을 보존하는, strain 과 독립인
    제3의 연속 손잡이. lambda=0 은 비극성 모체 Pbcm (대조군, Weyl 없어야 함).

(2) 벽 레지스트리  tau = 0 vs 0.5
    Pca2_1 의 두 glide 는 병진 성분에 둘 다 c 방향 1/2 을 가지므로, glide 는
    c 높이 h 의 벽을 h+1/2 의 벽으로 보낸다. 두 구조의 **총에너지 차이**가
    "자연이 한 레지스트리를 고르는가"를 결정한다 -- 이 갈래 전체가 여기 걸려 있다.

사용법:  python3 scripts/make_structures.py [--out structures]
"""

import argparse
import itertools
import os

import numpy as np


def read_poscar(path):
    L = open(path).read().splitlines()
    s = float(L[1].split()[0])
    cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    names = L[5].split()
    cnt = [int(x) for x in L[6].split()]
    assert L[7].strip()[0] in 'Dd', 'Direct 좌표만 지원'
    pos = np.array([[float(x) for x in L[8 + i].split()[:3]] for i in range(sum(cnt))])
    return cell, pos, names, cnt


def write_poscar(path, cell, pos, names, cnt, title):
    with open(path, 'w') as f:
        f.write(title + '\n   1.00000000000000\n')
        for v in cell:
            f.write('  %22.16f%22.16f%22.16f\n' % tuple(v))
        f.write('  ' + '  '.join(names) + '\n')
        f.write('  ' + '  '.join(str(c) for c in cnt) + '\n')
        f.write('Direct\n')
        for p in pos % 1.0:
            f.write('  %20.16f%20.16f%20.16f\n' % tuple(p))


def species(cnt):
    return np.concatenate([[i] * c for i, c in enumerate(cnt)]).astype(int)


def match(cell, z, pa, pb):
    from scipy.optimize import linear_sum_assignment
    idx = np.arange(len(pa))
    out = np.zeros_like(pb)
    worst = 0.0
    for Z in np.unique(z):
        ia = idx[z == Z]
        C = np.zeros((len(ia), len(ia)))
        for x, i in enumerate(ia):
            d = pb[ia] - pa[i]
            d -= np.round(d)
            C[x] = np.linalg.norm(d @ cell, axis=1)
        r, c = linear_sum_assignment(C)
        for x, y in zip(r, c):
            out[ia[x]] = pb[ia[y]]
            worst = max(worst, C[x, y])
    return out, worst


def min_dist(cell, pos):
    m = 9e9
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            d = pos[i] - pos[j]
            d -= np.round(d)
            m = min(m, np.linalg.norm(d @ cell))
    return m


def sg(cell, pos, z, symprec=1e-3):
    import spglib
    ds = spglib.get_symmetry_dataset((cell, pos % 1.0, z + 1), symprec=symprec)
    inv = any(np.array_equal(r, -np.eye(3, dtype=int)) for r in ds.rotations)
    return ds.international, len(ds.rotations), inv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='source/pristine_mirror/POSCAR')
    ap.add_argument('--out', default='structures')
    ap.add_argument('--polar', type=int, default=1, help='편극축 지표 (0-based)')
    ap.add_argument('--stack', type=int, default=2, help='벽 법선 지표 (0-based)')
    a = ap.parse_args()

    cell, pos, names, cnt = read_poscar(a.src)
    z = species(cnt)

    # --- 180도 짝을 만드는 거울 원점 찾기 -> 극성 모드 u, 비극성 모체 r0
    best = None
    for s in np.linspace(0, 1, 401, endpoint=False):
        q = pos.copy()
        q[:, a.polar] = (s - q[:, a.polar]) % 1.0
        pb, w = match(cell, z, pos, q)
        if best is None or w < best[1]:
            best = (s, w, pb)
    s, w, pb = best
    d = pb - pos
    d -= np.round(d)
    r0 = (pos + d / 2) % 1.0
    u = -d / 2
    print('거울 원점 s = %.4f, 최대 원자변위 %.4f A, 극성 모드 RMS %.4f A'
          % (s, w, np.sqrt(((u @ cell) ** 2).sum(1).mean())))

    d1 = os.path.join(a.out, 'polar_amplitude')
    d2 = os.path.join(a.out, 'wall_registry')
    os.makedirs(d1, exist_ok=True)
    os.makedirs(d2, exist_ok=True)

    print('\n[1] 극성 진폭 계열  (12원자, 셀 고정)')
    print('  lambda    x     파일                    공간군      ops  반전  최소거리')
    for lam in (1.00, 0.80, 0.60, 0.40, 0.00):
        q = (r0 + lam * u) % 1.0
        fn = os.path.join(d1, 'POSCAR_lambda%.2f' % lam)
        write_poscar(fn, cell, q, names, cnt,
                     'HfO2 polar amplitude lambda=%.2f (x=%.2f mirror fraction)' % (lam, (1 - lam) / 2))
        n, o, inv = sg(cell, q, z)
        print('   %.2f   %.3f   %-22s %-9s %3d  %-5s %.3f A'
              % (lam, (1 - lam) / 2, os.path.basename(fn), n, o, inv, min_dist(cell, q)))

    print('\n[2] 벽 레지스트리  (c 법선 180도 벽, 미이완)')
    print('  nA:nB  tau    원자수  파일                        공간군    ops  최소거리  주기')
    pB = pos.copy()
    pB[:, a.polar] = (s - pB[:, a.polar]) % 1.0
    for nA, nB in ((2, 1), (4, 2)):
        N = nA + nB
        c2 = cell.copy()
        c2[a.stack] = cell[a.stack] * N
        for tau in (0.0, 0.5):
            q = []
            for j in range(N):
                src = pos.copy() if j < nA else pB.copy()
                if j >= nA:
                    src[:, a.stack] = (src[:, a.stack] + tau) % 1.0
                src[:, a.stack] = (src[:, a.stack] + j) / N
                q.append(src)
            q = np.vstack(q)
            order = np.argsort(np.concatenate([z] * N), kind='stable')
            q = q[order]
            zz = np.concatenate([z] * N)[order]
            fn = os.path.join(d2, 'POSCAR_tau%.2f_%dto%d' % (tau, nA, nB))
            write_poscar(fn, c2, q, names, [c * N for c in cnt],
                         'HfO2 180deg domain superlattice %d:%d, wall registry tau=%.2f' % (nA, nB, tau))
            n, o, inv = sg(c2, q, zz)
            print('   %d:%d   %.2f   %4d   %-26s %-9s %3d  %.3f A   %.2f A'
                  % (nA, nB, tau, len(q), os.path.basename(fn), n, o,
                     min_dist(c2, q), np.linalg.norm(c2[a.stack])))

    print('\n생성 위치: %s/  (%s 는 대조군)' % (a.out, 'lambda1.00 = pristine, lambda0.00 = 비극성 모체'))


if __name__ == '__main__':
    main()
