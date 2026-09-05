#!/usr/bin/env python3
"""
bulkgap_scan.py — Simphony `BulkGap_cube_calc` 로 전역/국소 스캔을 돌리고 결과를
환산좌표·THz 로 정리한다. coarse 전역 -> 국소 refine 을 반복하는 용도.

전제: `patch/apply_gapcube_loto_fix.py` 가 적용된 Simphony.
     (미적용이면 gapshape3D 가 LO-TO 를 건너뛰어 엉뚱한 밴드를 본다.)

단위: GapCube.dat 의 Ev/Ec 는 THz 가 아니라 동역학 행렬 고윳값
      `w^2(THz^2) * 0.036749` (eV->Hartree) 다. 따라서 파일의 gap 은 dw 가 아니라
      d(w^2) 이고, 노드 근방에서
          gap_file ~= 0.073498 * w_bar(THz) * dw(THz)
      이 스크립트는 이를 되돌려 dw 를 THz 로 출력한다.

사용법:
    # 전역 coarse (기약 쐐기 전체)
    python3 scripts/bulkgap_scan.py --dir work/pp --src source/parent_pristine -N 81

    # 국소 refine (원점과 변 길이를 환산좌표로)
    python3 scripts/bulkgap_scan.py --dir work/pp --src source/parent_pristine \
        --origin 0.13 0.06 -0.01 --span 0.04 0.04 0.03 -N 61
"""

import argparse
import os
import subprocess

import numpy as np

HARTREE = 0.036749          # eV -> Hartree, gapshape3D 의 고윳값 배율


def recip(poscar):
    L = open(poscar).read().splitlines()
    s = float(L[1].split()[0])
    cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * s
    return cell, 2 * np.pi * np.linalg.inv(cell).T


def sym_indices(src, symprec=1e-4):
    """거울이 되는 지표와 편극축 지표 (0-based)."""
    import spglib
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=src + '/POSCAR', interface_mode='vasp')
    ds = spglib.get_symmetry_dataset((u.cell, u.scaled_positions, u.numbers), symprec=symprec)
    mir = set()
    for r in ds.rotations:
        if np.count_nonzero(r - np.diag(np.diag(r))):
            continue
        d = np.diag(r)
        if np.sum(d == -1) == 1:
            mir.add(int(np.where(d == -1)[0][0]))
    polar = [i for i in range(3) if i not in mir]
    return ds.international, sorted(mir), polar


def patch_pnin(path, origin, vecs, NN, thresh):
    L = open(path).read().splitlines()
    out, i = [], 0
    o, V = origin, vecs
    while i < len(L):
        x = L[i]
        t = x.strip()
        if t == 'KCUBE_BULK':
            out += ['KCUBE_BULK', ' %.8f %.8f %.8f' % tuple(o)] + \
                   [' %.8f %.8f %.8f' % tuple(v) for v in V]
            i += 5
            continue
        if t.startswith('BulkGap_cube_calc'):
            out.append('  BulkGap_cube_calc     = T'); i += 1; continue
        if t.startswith(('BulkBand_calc', 'WeylChirality_calc', 'SlabBand_calc',
                         'FindNodes_calc', 'Wanniercenter_calc')):
            out.append('  ' + t.split('=')[0].strip() + ' = F'); i += 1; continue
        if t.startswith('NumOccupied'):
            out.append('  NumOccupied = 17'); i += 1; continue
        if 'Nk1 = ' in x: out.append('  Nk1 = %d' % NN[0]); i += 1; continue
        if 'Nk2 = ' in x: out.append('  Nk2 = %d' % NN[1]); i += 1; continue
        if 'Nk3 = ' in x: out.append('  Nk3 = %d' % NN[2]); i += 1; continue
        if 'LOTO_method' in x: i += 1; continue          # 아래에서 강제로 다시 넣는다
        if 'Gap_threshold' in x:
            # phonopy 유래 hr.dat 이면 반드시 'phonopy'. 빠지면 기본값 'qe' 로 돌아가
            # Ewald 분할이 달라져 밴드가 0.03 THz 규모로 어긋난다 (실측).
            out.append("  LOTO_method = 'phonopy'")
            out.append('  Gap_threshold = %g' % thresh); i += 1; continue
        out.append(x); i += 1
    open(path, 'w').write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True, help='POSCAR/BORN/hr.dat/pn.in 이 있는 작업 폴더')
    ap.add_argument('--src', required=True, help='대칭 판정용 원본 폴더 (POSCAR)')
    ap.add_argument('--pn', default=None, help='pn.x 경로 (기본: <repo>/simphony/src/pn.x)')
    ap.add_argument('--origin', type=float, nargs=3, default=[0.0, 0.0, 0.0])
    ap.add_argument('--span', type=float, nargs=3, default=[0.5, 0.5, 0.5])
    ap.add_argument('--v1', type=float, nargs=3, default=None,
                    help='KCUBE 변 벡터를 직접 지정 (--span 대신). 평면을 훑을 때는 '
                         '얇은 벡터를 반드시 v3 에 둘 것 - v2 에 두면 Simphony 가 '
                         '엉뚱한 격자를 만든다 (검증됨).')
    ap.add_argument('--v2', type=float, nargs=3, default=None)
    ap.add_argument('--v3', type=float, nargs=3, default=None)
    ap.add_argument('-N', type=int, default=81)
    ap.add_argument('-N1', type=int, default=None, help='축별 격자수 (미지정시 -N)')
    ap.add_argument('-N2', type=int, default=None)
    ap.add_argument('-N3', type=int, default=None)
    ap.add_argument('--thresh', type=float, default=0.05, help='Gap_threshold (파일 단위)')
    ap.add_argument('--tol', type=float, default=None,
                    help='대칭면 판정 여유 (환산). 기본: 격자간격의 1.5배 '
                         '- 대칭면 바로 옆 격자점(nodal line/plane 꼬리)을 걸러낸다.')
    ap.add_argument('--no-run', action='store_true',
                    help='pn.x 를 다시 돌리지 않고 기존 GapCube.dat 만 재해석한다.')
    ap.add_argument('--top', type=int, default=12)
    a = ap.parse_args()

    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    pn = a.pn or os.path.join(root, 'simphony', 'src', 'pn.x')

    NN = [a.N1 or a.N, a.N2 or a.N, a.N3 or a.N]
    VEC = [a.v1, a.v2, a.v3]
    if any(v is None for v in VEC):
        VEC = [[a.span[0], 0, 0], [0, a.span[1], 0], [0, 0, a.span[2]]]
    VEC = [np.array(v, float) for v in VEC]
    name, mir, polar = sym_indices(a.src)
    print('%s  거울지표=%s  편극축=%s' % (name, [i + 1 for i in mir], [i + 1 for i in polar]))

    gc = os.path.join(a.dir, 'GapCube.dat')
    if not a.no_run:
        patch_pnin(os.path.join(a.dir, 'pn.in'), a.origin, VEC, NN, a.thresh)
        if os.path.exists(gc):
            os.remove(gc)
        subprocess.run([pn], cwd=a.dir, stdout=open(os.path.join(a.dir, 'run.log'), 'w'),
                       stderr=subprocess.STDOUT, check=False)
    if not os.path.exists(gc):
        print('GapCube.dat 이 없습니다. run.log 확인.')
        return

    d = np.loadtxt(gc, skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    _, B = recip(os.path.join(a.dir, 'POSCAR'))
    K = np.array([np.linalg.solve(B.T, r) for r in d[:, 0:3]])
    wbar = np.sqrt(np.maximum(d[:, 4], 0) / HARTREE)          # THz
    dw = d[:, 3] / (2 * HARTREE * np.maximum(wbar, 1e-12))    # THz

    step = max((np.linalg.norm(VEC[i]) / (NN[i] - 1.0)) for i in range(3) if NN[i] > 1)
    tol = a.tol if a.tol is not None else 1.5 * step

    def enforced(k):
        for i in mir:
            if abs(k[i]) < tol or abs(abs(k[i]) - 0.5) < tol:
                return True
        for i in polar:                                        # 편극축 nodal plane
            if abs(abs(k[i]) - 0.5) < tol:
                return True
        return False

    keep = np.array([not enforced(k) for k in K])
    print('격자 %dx%dx%d (최대간격 %.5f), 대칭면 여유 %.5f, 기록 %d개 -> 자유 %d개'
          % (NN[0], NN[1], NN[2], step, tol, len(d), keep.sum()))
    if keep.sum() == 0:
        print('  -> 후보 없음')
        return
    P, G, W = K[keep], dw[keep], wbar[keep]
    order = np.argsort(G)
    shown = []
    print('  최저 후보 (중복 클러스터 제외, dw 단위 THz)')
    for i in order:
        if any(np.linalg.norm(P[i] - s) < 0.02 for s in shown):
            continue
        shown.append(P[i])
        print('   k=(%+.5f, %+.5f, %+.5f)  dw=%.3e THz   f=%.4f THz (%.2f meV)'
              % (P[i][0], P[i][1], P[i][2], G[i], W[i], W[i] * 4.135667))
        if len(shown) >= a.top:
            break


if __name__ == '__main__':
    main()
