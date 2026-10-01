#!/usr/bin/env python3
"""
make_biaxial_vasp.py — biaxial strain 계열의 VASP 입력 일습을 만든다.

구조 하나당 3단계다.

  01_strain  변형 축을 고정한 셀 완화        (INCAR_strain)
  02_born    Born 유효전하 + eps_inf, DFPT   (INCAR_born)
  03_force   phonopy 유한변위 초격자 힘      (INCAR_force)

**축 규약을 가정하지 않는다.** 구조마다 어느 지표가 편극축인지가 다르므로
(CLAUDE.md 2.1), spglib 으로 직접 찾아 거울축 둘에만 변형을 걸고 편극축은
자유 relax 로 둔다. LATTICE_CONSTRAINTS 줄도 그에 맞춰 써 준다.

사용법
------
    python3 scripts/make_biaxial_vasp.py --ref source/pristine_mirror \
        --strains -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 --out vasp/biaxial/runs

    # 축 고정을 손으로 지정하고 싶을 때 (자동 판정을 덮어쓴다)
    python3 scripts/make_biaxial_vasp.py --ref source/pristine_mirror \
        --strains -1 --lattice-constraints ".TRUE. .FALSE. .TRUE."

POTCAR 는 만들지 않는다. 각 디렉터리에 직접 넣거나 --potcar 로 경로를 주면
심볼릭 링크를 걸어 준다.
"""

import argparse
import os
import shutil

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TPL = os.path.join(os.path.dirname(HERE), 'vasp', 'biaxial')

# 이 프로젝트 자체 데이터에서 뽑은 실효 푸아송비 (초기 추정값 전용).
# results/weyl_trend/SUMMARY.md 의 실측 strain 표에서
#   m1 -0.670, p1 -0.608, m2.5 -0.742, m3 -0.711  ->  평균 0.68
# 완화가 어차피 제 값을 찾으므로 이온 스텝 수만 줄이는 용도다.
POISSON = 0.68


def read_poscar(path):
    L = open(path).read().splitlines()
    scale = float(L[1].split()[0])
    cell = np.array([[float(x) for x in L[i].split()[:3]] for i in (2, 3, 4)]) * scale
    i = 5
    names = L[i].split(); i += 1
    counts = [int(x) for x in L[i].split()]; i += 1
    if L[i].strip()[0] in 'Ss':          # Selective dynamics
        i += 1
    mode = L[i].strip(); i += 1
    n = sum(counts)
    pos = np.array([[float(x) for x in L[i + j].split()[:3]] for j in range(n)])
    if mode[0] in 'Cc Kk':
        pos = pos @ np.linalg.inv(cell)
    return L[0], cell, names, counts, pos


def write_poscar(path, title, cell, names, counts, pos):
    with open(path, 'w') as f:
        f.write(title.rstrip() + '\n   1.00000000000000\n')
        for v in cell:
            f.write('  %22.16f %22.16f %22.16f\n' % tuple(v))
        f.write('  ' + '  '.join(names) + '\n')
        f.write('  ' + '  '.join(str(c) for c in counts) + '\n')
        f.write('Direct\n')
        for p in pos:
            f.write('  %20.16f %20.16f %20.16f\n' % tuple(p))


def axes(poscar):
    """(공간군, 편극축 0-based, 거울축 둘 0-based) 를 spglib 에서 뽑는다."""
    import spglib
    _, cell, names, counts, pos = read_poscar(poscar)
    numbers = []
    for k, c in enumerate(counts):
        numbers += [k + 1] * c
    ds = spglib.get_symmetry_dataset((cell, pos, numbers), symprec=1e-4)
    mirrors = set()
    for r in ds.rotations:
        if np.count_nonzero(r - np.diag(np.diag(r))):
            continue
        dg = np.diag(r)
        if np.sum(dg == -1) == 1:
            mirrors.add(int(np.where(dg == -1)[0][0]))
    polar = [i for i in range(3) if i not in mirrors]
    if len(polar) != 1:
        raise SystemExit('편극축이 %s 로 하나가 아니다. 구조를 확인할 것.' % polar)
    return ds.international, polar[0], sorted(mirrors)


def kpoints(path, mesh, comment):
    with open(path, 'w') as f:
        f.write('%s\n0\nGamma\n  %d %d %d\n  0 0 0\n' % (comment, *mesh))


def put_incar(src, dst, constraints=None):
    s = open(os.path.join(TPL, src)).read()
    if constraints is not None:
        out = []
        for line in s.splitlines():
            if line.strip().startswith('LATTICE_CONSTRAINTS'):
                line = 'LATTICE_CONSTRAINTS = %s' % constraints
            out.append(line)
        s = '\n'.join(out) + '\n'
    open(dst, 'w').write(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True, help='기준 구조 디렉터리 또는 POSCAR 경로')
    ap.add_argument('--strains', type=float, nargs='+', required=True,
                    help='biaxial strain, 퍼센트 단위. 예: -3 -2 -1 0 1')
    ap.add_argument('--out', default='vasp/biaxial/runs')
    ap.add_argument('--dim', default='2 2 2', help='phonopy 초격자 (기본 2 2 2)')
    ap.add_argument('--kmesh', type=int, nargs=3, default=[6, 6, 6], help='단위포 완화용')
    ap.add_argument('--kmesh-born', type=int, nargs=3, default=[8, 8, 8],
                    help='Born/eps_inf 용. 수렴이 느려 더 조밀하게')
    ap.add_argument('--poisson', type=float, default=POISSON,
                    help='편극축 초기 추정 수축비 (완화가 덮어쓴다)')
    ap.add_argument('--lattice-constraints', default=None,
                    help='자동 판정을 덮어쓴다. 예: ".TRUE. .FALSE. .TRUE."')
    ap.add_argument('--potcar', default=None, help='있으면 각 디렉터리에 링크')
    ap.add_argument('--keep-cell-noise', action='store_true',
                    help='격자벡터의 미세 비대각 잡음을 그대로 둔다 (기본은 0 으로 정리)')
    args = ap.parse_args()

    ref = args.ref
    if os.path.isdir(ref):
        ref = os.path.join(ref, 'POSCAR')
    sg, ip, mir = axes(ref)
    title, cell, names, counts, pos = read_poscar(ref)
    if not args.keep_cell_noise:
        tol = 1e-6 * np.linalg.norm(cell, axis=1)[:, None]
        n_zeroed = int((np.abs(cell) < tol).sum() - (cell == 0).sum())
        if n_zeroed:
            mx = np.abs(cell[np.abs(cell) < tol]).max()
            print('격자벡터의 미세 성분 %d개를 0 으로 정리했다 (최대 %.2e Å).'
                  '  그대로 두려면 --keep-cell-noise' % (n_zeroed, mx))
        cell = np.where(np.abs(cell) < tol, 0.0, cell)
    L = np.linalg.norm(cell, axis=1)

    auto = ['.FALSE.'] * 3
    auto[ip] = '.TRUE.'                      # 편극축만 자유
    auto = ' '.join(auto)
    lc = args.lattice_constraints or auto

    print('기준 구조 : %s   (%s)' % (ref, sg))
    print('  a1,a2,a3 = %.5f, %.5f, %.5f Å' % tuple(L))
    print('  편극축   = a%d  (%.5f Å)  -> 자유 relax' % (ip + 1, L[ip]))
    print('  거울축   = a%d, a%d           -> biaxial 변형 후 고정' % (mir[0] + 1, mir[1] + 1))
    print('  LATTICE_CONSTRAINTS = %s   (자동 판정: %s)' % (lc, auto))
    if args.lattice_constraints and lc != auto:
        print('  !! 경고: 지정값이 자동 판정과 다르다.')
        print('     VASP 에서 .TRUE. 는 "그 격자벡터를 완화" 다.')
        print('     위 지정대로면 a%s 가 고정되고 a%s 가 풀린다.'
              % (','.join(str(i + 1) for i, v in enumerate(lc.split()) if v.upper().startswith('.F')),
                 ','.join(str(i + 1) for i, v in enumerate(lc.split()) if v.upper().startswith('.T'))))
        print('     의도한 biaxial(거울축 고정 / 편극축 자유) 과 맞는지 확인할 것.')

    dim = [int(x) for x in args.dim.split()]
    ksup = [max(1, int(round(k / d))) for k, d in zip(args.kmesh, dim)]
    print('  k메시 : 완화 %dx%dx%d / Born %dx%dx%d / 초격자(%s) %dx%dx%d'
          % (*args.kmesh, *args.kmesh_born, args.dim, *ksup))

    for e in args.strains:
        tag = 'eps_%s%05.2f' % ('m' if e < 0 else 'p', abs(e))
        root = os.path.join(args.out, tag)
        c = cell.copy()
        for i in mir:
            c[i] *= (1 + e / 100.0)                       # 변형 축
        c[ip] *= (1 - args.poisson * e / 100.0)           # 편극축 초기 추정
        d1 = os.path.join(root, '01_strain')
        d2 = os.path.join(root, '02_born')
        d3 = os.path.join(root, '03_force')
        for d in (d1, d2, d3):
            os.makedirs(d, exist_ok=True)

        write_poscar(os.path.join(d1, 'POSCAR'),
                     '%s  biaxial %+.2f%%  (a%d,a%d strained+fixed; a%d free)'
                     % (' '.join(title.split()[:2]), e, mir[0] + 1, mir[1] + 1, ip + 1),
                     c, names, counts, pos)
        put_incar('INCAR_strain', os.path.join(d1, 'INCAR'), constraints=lc)
        put_incar('INCAR_born', os.path.join(d2, 'INCAR'))
        put_incar('INCAR_force', os.path.join(d3, 'INCAR'))
        kpoints(os.path.join(d1, 'KPOINTS'), args.kmesh, 'relax %+.2f%%' % e)
        kpoints(os.path.join(d2, 'KPOINTS'), args.kmesh_born, 'born %+.2f%%' % e)
        kpoints(os.path.join(d3, 'KPOINTS'), ksup, 'force supercell %s' % args.dim)
        if args.potcar:
            for d in (d1, d2, d3):
                lnk = os.path.join(d, 'POTCAR')
                if not os.path.exists(lnk):
                    os.symlink(os.path.abspath(args.potcar), lnk)

        open(os.path.join(root, 'run.sh'), 'w').write(RUNSH % dict(
            tag=tag, dim=args.dim, eps=e))
        os.chmod(os.path.join(root, 'run.sh'), 0o755)
        print('  %-12s a=(%.5f, %.5f, %.5f) 초기추정' % (tag, *np.linalg.norm(c, axis=1)))

    print('\n각 디렉터리의 run.sh 가 세 단계를 순서대로 돈다. POTCAR 는 직접 넣을 것.')


RUNSH = r'''#!/bin/bash
# %(tag)s : biaxial %(eps)+.2f%%   — 세 단계를 순서대로
set -e
VASP=${VASP:-vasp_std}
MPI=${MPI:-"mpirun -np ${NPROC:-16}"}
ROOT=$(cd "$(dirname "$0")" && pwd)

# --- 1. 축 고정 완화 (Pulay 응력 때문에 2회 반복) ---------------------
cd "$ROOT/01_strain"
for i in 1 2; do
  $MPI $VASP
  cp CONTCAR POSCAR
  cp OUTCAR OUTCAR.relax$i
done
echo "완화 완료. 최종 격자:"; sed -n '3,5p' POSCAR

# --- 2. Born 유효전하 (완화된 단위포) --------------------------------
cd "$ROOT/02_born"
cp ../01_strain/CONTCAR POSCAR
$MPI $VASP
phonopy-vasp-born > BORN
echo "BORN 생성:"; cat BORN

# --- 3. 유한변위 힘 ---------------------------------------------------
cd "$ROOT/03_force"
cp ../01_strain/CONTCAR POSCAR
phonopy -d --dim="%(dim)s" -c POSCAR
for f in POSCAR-*; do
  d="disp-${f#POSCAR-}"
  mkdir -p "$d"
  cp "$f" "$d/POSCAR"
  cp INCAR KPOINTS POTCAR "$d/" 2>/dev/null || true
  ( cd "$d" && $MPI $VASP )
done
phonopy -f disp-*/vasprun.xml
echo "FORCE_SETS 생성 완료"

# --- 결과 모으기 ------------------------------------------------------
mkdir -p "$ROOT/out"
cp "$ROOT/01_strain/CONTCAR" "$ROOT/out/POSCAR"
cp "$ROOT/02_born/BORN"      "$ROOT/out/BORN"
cp "$ROOT/03_force/FORCE_SETS" "$ROOT/out/FORCE_SETS"
echo "out/ 에 POSCAR, BORN, FORCE_SETS 가 모였다. weyl_scan.py 에 바로 쓸 수 있다."
'''

if __name__ == '__main__':
    main()
