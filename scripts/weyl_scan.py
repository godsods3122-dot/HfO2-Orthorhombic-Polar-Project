#!/usr/bin/env python3
"""
weyl_scan.py — phonopy를 참값으로 써서 Weyl point / nodal line을 찾는다.

이 대화에서 인라인으로 돌렸던 탐색 과정을 하나로 합친 것. 패치된 Simphony와
전 경로 3.4e-7 THz로 일치하므로(README 참조), 노드 위치 탐색은 이쪽이 훨씬 빠르다.
Simphony는 chirality 계산에만 쓰면 된다.

핵심 원리
---------
1) 대칭 조작을 spglib에서 직접 뽑는다. 구조마다 축 규약이 달라서
   "어느 지표가 거울면인가"가 바뀐다. 이걸 가정하면 반드시 틀린다.
   - 거울 조작(R의 대각이 하나만 -1)의 불변면 k_i = 0, 0.5 위에서는
     chirality가 0으로 강제된다. Weyl은 여기 살 수 없다.
   - 편극축(2_1 screw 방향) 지표에는 거울이 없으므로 Weyl이 살 수 있다.
2) 기약 쐐기만 훑는다. 점군이 (±k_a, ±k_b, k_c) 꼴이면 [0,0.5]^3로 충분하다.
3) 격자에서 국소 최소를 뽑고 Nelder-Mead로 정밀화한 뒤,
   거울면 위로 수렴하는 것은 버린다.

사용법
------
    python weyl_scan.py --dir . --band 17 --mode sym      # 대칭 조작만 출력
    python weyl_scan.py --dir . --band 17 --mode wedge -N 41
    python weyl_scan.py --dir . --band 17 --mode plane --fixed-index 2 --fixed-value 0.0 -N 161
    python weyl_scan.py --dir . --band 17 --mode refine --seed 0.116 0.0 0.005

--band 17 이면 17번과 18번 밴드 사이의 gap을 본다 (0-based 16,17).
--dir 에는 POSCAR, FORCE_SETS, BORN 이 있어야 한다.
"""

import argparse
import itertools
import sys

import numpy as np


def get_ph(d, dim=(2, 2, 2)):
    from phonopy import Phonopy
    from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    ph = Phonopy(u, supercell_matrix=np.diag(dim), primitive_matrix='P')
    ph.dataset = parse_FORCE_SETS(filename=d + '/FORCE_SETS')
    ph.produce_force_constants(calculate_full_force_constants=True)
    ph.symmetrize_force_constants_by_space_group()
    ph.symmetrize_force_constants(level=3)
    P = parse_BORN(ph.primitive, filename=d + '/BORN')
    P['factor'] = 14.399652          # eV*Angstrom, phonopy 기본값
    P['method'] = 'gonze'
    ph.nac_params = P
    return ph


def symmetry_info(d):
    """거울면이 되는 지표와 편극축 지표를 반환한다."""
    import spglib
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    ds = spglib.get_symmetry_dataset((u.cell, u.scaled_positions, u.numbers), symprec=1e-4)
    mirrors = set()
    for r in ds.rotations:
        if np.count_nonzero(r - np.diag(np.diag(r))):
            continue                      # 대각이 아닌 조작은 여기서 다루지 않음
        dg = np.diag(r)
        if np.sum(dg == -1) == 1:         # 지표 하나만 뒤집으면 거울
            mirrors.add(int(np.where(dg == -1)[0][0]))
    polar = [i for i in range(3) if i not in mirrors]
    return ds.international, sorted(mirrors), polar, ds.rotations, ds.translations


def gap_at(ph, qs, band):
    ph.run_qpoints([list(q) for q in qs])
    f = ph.qpoints.frequencies
    return f[:, band] - f[:, band - 1]


def refine(ph, seed, band):
    from scipy.optimize import minimize
    def g(q):
        return float(gap_at(ph, [q], band)[0])
    r = minimize(g, list(seed), method='Nelder-Mead',
                 options=dict(xatol=1e-9, fatol=1e-14, maxiter=8000, maxfev=8000))
    return r.x, r.fun


def on_mirror(k, mirrors, tol=1e-4):
    for i in mirrors:
        if abs(k[i]) < tol or abs(abs(k[i]) - 0.5) < tol:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='.')
    ap.add_argument('--band', type=int, default=17, help='이 밴드와 다음 밴드 사이 gap')
    ap.add_argument('--dim', type=int, nargs=3, default=[2, 2, 2])
    ap.add_argument('--mode', choices=['sym', 'wedge', 'plane', 'refine'], default='wedge')
    ap.add_argument('-N', type=int, default=41)
    ap.add_argument('--thresh', type=float, default=0.02, help='국소최소 후보 gap 상한 (THz)')
    ap.add_argument('--fixed-index', type=int, default=2)
    ap.add_argument('--fixed-value', type=float, default=0.0)
    ap.add_argument('--seed', type=float, nargs=3)
    args = ap.parse_args()

    name, mirrors, polar, R, T = symmetry_info(args.dir)
    print("공간군: %s" % name)
    print("거울면이 되는 지표 (여기서는 chirality=0 강제): %s" % [i + 1 for i in mirrors])
    print("편극축 지표 (Weyl 가능)                      : %s" % [i + 1 for i in polar])
    for i, (r, t) in enumerate(zip(R, T)):
        print("   op%d  R대각=%s  t=%s" % (i, np.diag(r), np.round(t, 3)))
    if args.mode == 'sym':
        return

    ph = get_ph(args.dir, tuple(args.dim))
    b = args.band

    if args.mode == 'refine':
        if not args.seed:
            sys.exit("--seed k1 k2 k3 를 주세요.")
        x, g = refine(ph, args.seed, b)
        tag = "거울면 위 (chirality 0 강제)" if on_mirror(x, mirrors) else "일반 위치 (Weyl 후보)"
        print("\nk=(%.7f, %.7f, %.7f)  gap=%.3e THz   -> %s" % (x[0], x[1], x[2], g, tag))
        return

    if args.mode == 'plane':
        fi, fv, N = args.fixed_index, args.fixed_value, args.N
        ax = [i for i in range(3) if i != fi]
        us = np.linspace(0.0, 0.5, N)
        G = np.zeros((N, N))
        for i, a in enumerate(us):
            qs = []
            for c in us:
                q = [0.0, 0.0, 0.0]
                q[fi] = fv; q[ax[0]] = a; q[ax[1]] = c
                qs.append(q)
            G[i] = gap_at(ph, qs, b)
        print("\nk%d=%.3f 평면 (간격 %.4f)" % (fi + 1, fv, us[1] - us[0]))
        print("  전체 최소 gap        = %.4e" % G.min())
        inner = G[3:-3, 3:-3]
        print("  경계에서 떨어진 최소 = %.4e" % inner.min())
        cnt = 0
        for i in range(1, N - 1):
            for j in range(1, N - 1):
                if G[i, j] < args.thresh and G[i, j] == G[i-1:i+2, j-1:j+2].min():
                    k = [0.0, 0.0, 0.0]; k[fi] = fv; k[ax[0]] = us[i]; k[ax[1]] = us[j]
                    if not on_mirror(k, mirrors):
                        print("   후보 k=(%.4f,%.4f,%.4f) gap=%.4e" % (k[0], k[1], k[2], G[i, j]))
                        cnt += 1
        print("  거울면 밖 국소최소: %d개" % cnt)
        return

    # wedge
    N = args.N
    ks = np.linspace(0.0, 0.5, N)
    G = np.zeros((N, N, N), dtype=np.float32)
    for i, a in enumerate(ks):
        qs = [[a, x, y] for x in ks for y in ks]
        G[i] = gap_at(ph, qs, b).reshape(N, N)
        if i % 10 == 0:
            print("  k1 %d/%d" % (i, N), flush=True)
    print("\n쐐기 최소 gap = %.4e" % G.min())

    cands = []
    for i in range(N):
        for j in range(N):
            for k in range(N):
                v = G[i, j, k]
                if v >= args.thresh:
                    continue
                sl = G[max(0, i-1):i+2, max(0, j-1):j+2, max(0, k-1):k+2]
                if v == sl.min():
                    cands.append((float(v), ks[i], ks[j], ks[k]))
    cands.sort()
    print("gap<%.3f 국소최소 %d개" % (args.thresh, len(cands)))

    weyl, trivial = [], 0
    for v, a, bb, c in cands:
        if on_mirror([a, bb, c], mirrors):
            trivial += 1
            continue
        x, g = refine(ph, [a, bb, c], b)
        if on_mirror(x, mirrors):
            trivial += 1
        else:
            weyl.append((g, x))
    print("  거울면으로 수렴/거울면 위 : %d개  (chirality 0 강제)" % trivial)
    print("  일반 위치에 남은 것       : %d개" % len(weyl))
    for g, x in sorted(weyl)[:20]:
        print("   k=(%.7f,%.7f,%.7f)  gap=%.3e" % (x[0], x[1], x[2], g))
    if weyl:
        np.savetxt('weyl_candidates.dat', np.array([list(x) + [g] for g, x in weyl]),
                   header='k1 k2 k3 gap', fmt='%18.10f')
        print("\nweyl_candidates.dat 저장. Simphony WEYL_CHIRALITY 카드에 넣어 확인하세요.")


if __name__ == '__main__':
    main()
