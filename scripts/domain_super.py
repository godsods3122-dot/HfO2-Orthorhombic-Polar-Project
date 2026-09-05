#!/usr/bin/env python3
"""
domain_super.py — 180° 강유전 도메인(= parent/mirror 쌍) 배열의 공간군을 계산한다.

동기
----
조성·strain 을 전혀 건드리지 않고, **각 도메인 안은 100% 원래 Pca2₁ HfO₂ 그대로**
두면서 초격자 수준에서만 대칭을 낮추는 길. 도메인 벽은 HfO₂ 의 고유 미세구조이고
전기장으로 쓰고 지울 수 있으므로 "물성을 바꾸지 않으면서 유동적으로 조절"이라는
요구를 그대로 만족한다.

핵심 결과 (results/domain_route/SUMMARY.md 참조)
-----------------------------------------------
  * 도메인 폭이 같으면(미분극/보상 상태) 초격자가 **반전대칭**을 얻는다
    → 시간역전과 함께 Berry curvature 를 전 BZ 에서 0 으로 만들어 **Weyl 이 전멸**.
  * 도메인 폭이 다르면(부분 분극) 반전이 없고, **변조 방향이 결정한다**:
      - 변조 ∥ a  →  Pca2₁ 그대로 (아무것도 안 깨짐, folding 만)
      - 변조 ∥ c  →  **P2₁** (거울 0개 — 표준 고대칭 경로가 열린다)
    이유: Pca2₁ 의 두 glide 는 병진 성분에 **둘 다 c 방향 ½** 을 갖는다
    (op2 t=(0,½,½), op3 t=(0,0,½)). 그래서 c 방향 변조는 둘 다 깨고,
    a 방향 변조는 둘 다 안 깬다.
  * 벽면이 편극축을 포함하므로(법선 ⊥ 편극축) 전기적으로 중성 = 물리적으로 선호되는
    180° 벽이다. 편극축 방향 변조(head-to-head)는 최소원자거리 0.09 Å 로 비물리적.

주의: 축 규약은 구조마다 다르다. 반드시 `weyl_scan.py --mode sym` 으로 먼저 확인할 것
(이 스크립트의 --polar 기본값은 pristine_mirror 기준).

사용법
------
    python3 scripts/domain_super.py --dir source/pristine_mirror
    python3 scripts/domain_super.py --dir source/m1_mirror --polar 1
"""

import argparse

import numpy as np


def load(d):
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    return np.array(u.cell), np.array(u.scaled_positions), np.array(u.numbers)


def flip(p, polar):
    """180° 도메인 짝 = 편극축 분율좌표 뒤집기 (레포에서 1.8e-5 로 검증된 관계)."""
    q = p.copy()
    q[:, polar] = (-q[:, polar]) % 1.0
    return q


def min_dist(L, q):
    m = 9e9
    for i in range(len(q)):
        for j in range(i + 1, len(q)):
            d = q[i] - q[j]
            d -= np.round(d)
            m = min(m, np.linalg.norm(d @ L))
    return m


def stack(L, p, pB, n, direction, nA, nB, symprec=1e-3):
    """도메인 A 를 nA 셀, B 를 nB 셀 쌓은 초격자의 spglib dataset 을 돌려준다."""
    import spglib
    N = nA + nB
    L2 = L.copy()
    L2[direction] = L[direction] * N
    q, num = [], []
    for j in range(N):
        src = p if j < nA else pB
        c = src.copy()
        c[:, direction] = (src[:, direction] + j) / N
        q.append(c)
        num.append(n)
    q = np.vstack(q)
    num = np.concatenate(num)
    return spglib.get_symmetry_dataset((L2, q, num), symprec=symprec), min_dist(L2, q)


def describe(ds):
    R = ds.rotations
    inv = any(np.array_equal(r, -np.eye(3, dtype=int)) for r in R)
    mir = sorted({int(np.where(np.diag(r) == -1)[0][0]) for r in R
                  if np.array_equal(np.diag(np.diag(r)), r) and sum(np.diag(r) == -1) == 1})
    return ds.international, len(R), inv, mir


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='source/pristine_mirror')
    ap.add_argument('--polar', type=int, default=1, help='편극축 지표 (0-based)')
    ap.add_argument('--symprec', type=float, default=1e-3)
    a = ap.parse_args()

    L, p, n = load(a.dir)
    pB = flip(p, a.polar)
    ax = 'abc'

    import spglib
    ds0 = spglib.get_symmetry_dataset((L, p, n), symprec=a.symprec)
    print('단일 도메인 (완전 분극): %s   거울지표=%s   net P = +1.000 P0'
          % (ds0.international, describe(ds0)[3]))
    print('편극축 지표 = %d (%s)\n' % (a.polar, ax[a.polar]))

    hdr = ' 변조   nA:nB   최소원자거리   공간군      ops  반전  거울지표   net P'
    for label, ratios in [('대칭 배열 (미분극 / 완전 보상)', [(1, 1), (2, 2), (4, 4)]),
                          ('비대칭 배열 (부분 분극)', [(2, 1), (3, 1), (5, 3), (7, 5), (8, 3)])]:
        print('== %s ==' % label)
        print(hdr)
        for D in range(3):
            for nA, nB in ratios:
                ds, md = stack(L, p, pB, n, D, nA, nB, a.symprec)
                sg, no, inv, mir = describe(ds)
                flag = '   <-- 비물리적 벽 (head-to-head)' if md < 1.0 else ''
                print('   %s   %d:%d    %6.3f A    %-9s %3d  %-5s %-9s %+0.3f%s'
                      % (ax[D], nA, nB, md, sg, no, inv, str(mir), (nA - nB) / (nA + nB), flag))
        print()


if __name__ == '__main__':
    main()
