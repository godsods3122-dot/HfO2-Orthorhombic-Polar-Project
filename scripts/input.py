#!/usr/bin/env python3
"""
Generate 'pn.in' for Phonnier from Phonopy (POSCAR + optional BORN).

Usage:
    python input.py --poscar POSCAR --born BORN --hr wannier90_hr.dat

phonopy V3 API 점검 결과 (2026-07 기준 공식 문서/phono3py API 문서와 대조):
  - parse_BORN(ph.primitive, filename=...) 형태는 현재 API와 일치 (변경 없음)
  - primitive.numbers / .symbols / .masses / .scaled_positions / .cell 는
    PhonopyAtoms의 표준 속성으로 현재도 유효 (변경 없음)
  - primitive_matrix='auto' 문자열 지정 자체는 V3에서도 유효 (변경 없음).
    단, phonopy 4.x부터는 이 값이 "기본값"으로 바뀌고 CLI가 분리되는 등
    breaking change가 있으므로, 4.x로 넘어가면 다시 점검이 필요함.

**실질적으로 고친 부분은 API가 아니라 다음의 일관성 문제입니다**:
  이전에 만든 phonopy2TBDAT.py(FORCE_SETS+POSCAR -> wannier90_hr.dat)는
  R벡터 정수화 버그 때문에 primitive_matrix 기본값을 'identity'로 씁니다.
  이 스크립트가 여전히 'auto'를 쓰면, hr.dat 파일과 pn.in의 LATTICE/ATOMS/
  KPATH_BULK/WEYL_CHIRALITY가 서로 다른(회전/재정렬된) primitive cell
  기준으로 쓰이게 되어 Phonnier 계산 전체가 어긋납니다. 또한 BORN 파일도
  특정 primitive cell을 기준으로 대칭적으로 독립인 원자만 나열되어 있으므로,
  parse_BORN에 넘기는 primitive가 BORN을 만들 때 쓴 것과 다르면 파싱 자체가
  조용히 어긋날 수 있습니다. 그래서 phonopy2TBDAT.py와 동일한
  --primitive-matrix 옵션(기본값 identity)을 추가해 두 스크립트가 항상 같은
  좌표계를 쓰도록 맞췄습니다. hr.dat을 만들 때 -p 옵션을 바꿨다면, 이
  스크립트도 반드시 동일한 값으로 실행하세요.
"""

import argparse
from datetime import datetime

import numpy as np
import seekpath

from phonopy import Phonopy
from phonopy.file_IO import parse_BORN
from phonopy.interface.calculator import read_crystal_structure


_BRAVAIS_KEYWORDS = {"p": "P", "f": "F", "i": "I", "a": "A", "c": "C", "r": "R"}


def parse_matrix_arg(s):
    nums = [float(x) for x in s.replace(",", " ").split()]
    if len(nums) == 3:
        return np.diag(nums)
    if len(nums) == 9:
        return np.array(nums).reshape(3, 3)
    raise argparse.ArgumentTypeError(
        "행렬 인자는 숫자 3개(대각) 또는 9개(3x3, row-major)여야 합니다: '{}'".format(s)
    )


def parse_primitive_matrix_arg(s):
    key = s.strip().lower()
    if key in ("identity", "none", "unit", "unitcell"):
        # phonopy2TBDAT.py의 기본값과 동일: POSCAR 자체를 primitive cell로
        # 그대로 사용한다. phonopy v4부터 None은 새 기본값 'auto'로 풀리므로
        # (PrimitiveMatrixAutoDefaultWarning) None이 아니라 명시적으로 'P'
        # (Bravais: 센터링 없음 = identity)를 넘겨 v3/v4 모두에서 identity를
        # 보장한다.
        return "P"
    if key == "auto":
        return "auto"
    if key in _BRAVAIS_KEYWORDS:
        return _BRAVAIS_KEYWORDS[key]
    return parse_matrix_arg(s)


def main():
    parser = argparse.ArgumentParser(description="Generate pn.in from Phonopy structure")
    parser.add_argument("--poscar", type=str, default="POSCAR", help="Path to unit cell POSCAR")
    parser.add_argument("--born", type=str, default=None, help="Path to BORN file for LOTO correction (optional)")
    parser.add_argument("--hr", type=str, default="wannier90_hr.dat", help="Name of the tight-binding HR file")
    parser.add_argument("--dim", type=str, default="1 1 1", help="Supercell dimensions for Phonopy init (doesn't affect primitive)")
    parser.add_argument(
        "-p", "--primitive-matrix", default="identity", type=parse_primitive_matrix_arg,
        help="primitive matrix. 'identity'(POSCAR 자체가 원시셀, 기본값. "
             "phonopy2TBDAT.py와 반드시 동일하게 맞출 것) / 'auto' / "
             "'P','F','I','A','C','R' / '3 또는 9개 숫자'"
    )
    parser.add_argument(
        "--transpose-born", dest="transpose_born", action="store_true", default=False,
        help="LOTO_BC에 쓰기 전에 각 원자의 Born charge 3x3 텐서를 전치한다 "
             "(기본값: 꺼짐). 2026-08-29 측정으로 비전치가 옳음이 확정됐다 - parent_m1 데이터셋에서 비전치는 전 36밴드에서 phonopy와 2.6e-7 THz 일치, 전치는 1.8e-1 THz 어긋나고 Weyl 자리에서 gap이 닫히지도 않는다. 아래 옛 근거는 기록용이다: Simphony 공식 Al2ZnTe4 예제의 LOTO_BC를 phonopy "
             "원본 Born charge와 원자별로 대조한 결과, Simphony는 전치된 "
             "텐서를 기대하는 것으로 확인됨."
    )
    parser.add_argument(
        "--no-transpose-born", dest="transpose_born", action="store_false",
        help="전치를 끄고 phonopy가 반환하는 원본 축 순서 그대로 쓴다."
    )
    args = parser.parse_args()

    print(f"Loading structure from {args.poscar}...")
    unitcell, _ = read_crystal_structure(args.poscar, interface_mode='vasp')

    dim_list = list(map(int, args.dim.split()))
    supercell_matrix = np.diag(dim_list) if len(dim_list) == 3 else np.array(dim_list).reshape(3, 3)

    ph = Phonopy(unitcell, supercell_matrix=supercell_matrix, primitive_matrix=args.primitive_matrix)
    primitive = ph.primitive

    cell = primitive.cell
    positions = primitive.scaled_positions
    numbers = primitive.numbers
    masses = primitive.masses
    labels = primitive.symbols
    num_atom = len(positions)

    isloto = "F"
    dt = None
    bc = None
    if args.born:
        try:
            print(f"Reading Born effective charges from {args.born}...")
            nac_params = parse_BORN(ph.primitive, filename=args.born)
            if nac_params:
                dt = nac_params['dielectric']
                bc = nac_params['born']
                if args.transpose_born:
                    print("--transpose-born: 각 원자의 Born charge 텐서를 전치합니다.")
                    bc = np.transpose(bc, axes=(0, 2, 1))

                # --- 검증 1: 배열 크기가 원시셀 원자 수와 일치하는가 ---
                # BORN 파일은 보통 phonopy-vasp-born 같은 별도 도구가 "그
                # 도구 자신의" 대칭 판단으로 독립 원자를 추려 만든다. 지금
                # ph.primitive의 대칭 분석(symprec, 구조 미세 차이 등)이
                # 그것과 다르면 parse_BORN이 원자 수를 잘못 확장할 수 있다.
                if bc.shape[0] != num_atom:
                    raise ValueError(
                        "BORN 파일에서 확장된 Born charge 원자 수({})가 "
                        "primitive cell의 원자 수({})와 다릅니다. BORN 파일을 "
                        "만들 때 쓴 구조/대칭 기준이 지금 --poscar/-p 설정과 "
                        "다른지 확인하세요.".format(bc.shape[0], num_atom)
                    )

                # --- 검증 2: 음향 총합 규칙 (acoustic sum rule) ---
                # 전체 Born charge의 합은 이상적으로 0(전하 중성)이어야 한다.
                # DFT BEC 계산이 완전히 대칭적이지 않거나 대칭 확장이 어긋나면
                # 이 합이 눈에 띄게 0에서 벗어난다 - 조용한 오적용을 잡아내는
                # 가장 간단한 체크.
                born_sum = np.sum(bc, axis=0)
                born_sum_norm = np.linalg.norm(born_sum)
                print(f"Born charge 합(음향 총합 규칙, 0에 가까워야 함): "
                      f"norm={born_sum_norm:.6f}")
                print(f"{born_sum}")
                if born_sum_norm > 0.1:
                    print(
                        "경고: Born charge 합이 0에서 크게 벗어났습니다. "
                        "독립 원자 -> 전체 원자 확장이 잘못됐을 가능성이 있으니 "
                        "BORN 파일과 POSCAR/primitive-matrix 조합을 다시 "
                        "확인하세요. 그래도 LOTO_correction은 켠 채로 "
                        "진행합니다 - 계산 전에 직접 판단하세요."
                    )

                isloto = "T"
                print("LOTO parameters successfully loaded.")
        except Exception as e:
            print(f"Failed to load BORN file: {e}. Proceeding without LOTO.")
            isloto = "F"
    else:
        print("--born 미지정: LOTO_correction = F 로 진행합니다 "
              "(hr.dat이 short-range 힘상수만 담고 있다는 전제와 일치).")

    structure = (cell, positions, numbers)
    path = seekpath.get_explicit_k_path_orig_cell(structure, with_time_reversal=True, reference_distance=0.05)

    qpath = path['explicit_kpoints_rel']
    segments = path['explicit_segments']
    qlabels = path['explicit_kpoints_labels']

    for i in range(len(qlabels)):
        if qlabels[i] == "GAMMA":
            qlabels[i] = "G"

    outfile = "pn.in"
    print(f"Writing {outfile}...")
    with open(outfile, "w") as f:
        f.write(f"!> Input file written automatically on {datetime.now()} by Phonopy-to-Phonnier script\n")

        f.write("&TB_FILE\n"
                f"  Hrfile = '{args.hr}'\n"
                f"  Package = 'Phonopy'\n"
                "/\n\n")

        f.write("&CONTROL\n"
                "!> LO-TO correction (if True, long range part of Dynamical matrix must be removed from phononTB_hr.dat)\n"
                f"  LOTO_correction       = {isloto}\n\n"
                "!> bulk band structure calculation flags\n"
                "  BulkBand_calc         = F\n"
                "  BulkGap_cube_calc     = F\n"
                "  BulkGap_plane_calc    = F\n\n"
                "!> slab band structure calculation flags\n"
                "  SlabBand_calc         = F\n"
                "  SlabSS_calc           = F\n"
                "  SlabArc_calc          = F\n\n"
                "!> wire band structure calculation flags\n"
                "  WireBand_calc         = F\n\n"
                "!> DOS calculation flags\n"
                "  Dos_calc              = F\n\n"
                "!> Topological quantities calculation flags\n"
                "  FindNodes_calc        = F\n"
                "  BerryPhase_calc       = F\n"
                "  BerryCurvature_calc   = F\n"
                "  Chern_3D_calc         = F\n"
                "  Wanniercenter_calc    = F\n"
                "  WeylChirality_calc    = F\n"
                "/\n\n")

        f.write("&SYSTEM \n"
                "  NSLAB = 10\n"
                "  NumOccupied = 1      !>Automatically selects from band 1 up to band NumOccupied. Overwriten by SELECTED_OCCUPIED_BANDS card\n"
                "/\n\n")

        f.write("&PARAMETERS \n"
                "  Eta_Arc = 0.005    !>infinite small value, like brodening. Too small a value leads to singularities in the Green Function \n"
                "  E_arc = 1      !>energy for calculate Fermi Arc, THz \n"
                "  OmegaNum = 100     !>omega number \n"
                "  OmegaMin =  0.5      !>energy interval in unit of THz \n"
                "  OmegaMax =  1.5     !>energy interval in unit of THz \n"
                "  Nk1 = 50           !>number k points \n"
                "  Nk2 = 50           !>number k points \n"
                "  Nk3 = 50          !>number k points \n"
                "  NP = 5             !>number of principle layers \n"
                "  Gap_threshold = 0.00001 !>threshold for GapCube output \n"
                "/\n\n")

        f.write("LATTICE\n"
                "Angstrom\n"
                f"{cell[0][0]:20.10f}    {cell[0][1]:20.10f}    {cell[0][2]:20.10f}\n"
                f"{cell[1][0]:20.10f}    {cell[1][1]:20.10f}    {cell[1][2]:20.10f}\n"
                f"{cell[2][0]:20.10f}    {cell[2][1]:20.10f}    {cell[2][2]:20.10f}\n\n")

        f.write("ATOMS\n"
                f"{num_atom}\n"
                "Direct\n")
        for i in range(num_atom):
            f.write(f"{labels[i]:4s}    {masses[i]:10.5f}    {positions[i][0]:20.10f}    {positions[i][1]:20.10f}    {positions[i][2]:20.10f}\n")
        f.write("\n")

        f.write("KPATH_BULK\n"
                f"{len(segments)}\n")

        for i, s in enumerate(segments):
            i1, i2 = s[0], s[1] - 1
            f.write(f"{qlabels[i1]:4s}    {qpath[i1][0]:20.8f}    {qpath[i1][1]:20.8f}    {qpath[i1][2]:20.8f}    "
                    f"{qlabels[i2]:4s}    {qpath[i2][0]:20.8f}    {qpath[i2][1]:20.8f}    {qpath[i2][2]:20.8f}\n")
        f.write("\n")

        f.write("SELECTED_OCCUPIED_BANDS\n"
                "1-3 !>Selects bands from 1 to 3\n\n")

        # 주의: 아래 SURFACE/WEYL_CHIRALITY 값은 예시 placeholder입니다.
        # 실제 계산된 Weyl point 좌표/표면 방향으로 반드시 교체하세요.
        f.write("SURFACE\n"
                "1 0 0\n"
                "0 1 0\n\n")

        f.write("KPATH_SLAB\n"
                "3\n"
                "G 0.0 0.0 X 0.5 0.0\n"
                "X 0.5 0.0 M 0.5 0.5\n"
                "M 0.5 0.5 G 0.0 0.0\n\n")

        f.write("KPLANE_SLAB\n"
                "-0.5 -0.5\n"
                "1.0 0.0\n"
                "0.0 1.0\n\n")

        f.write("KPLANE_BULK\n"
                "-0.00 -0.00 -0.00\n"
                "1.00 0.00 0.00\n"
                "0.00 1.00 0.00\n\n")

        f.write("KCUBE_BULK\n"
                "-0.00 -0.00 -0.00\n"
                "1.00 0.00 0.00\n"
                "0.00 1.00 0.00\n"
                "0.00 0.00 1.00\n\n")

        f.write("WEYL_CHIRALITY\n"
                "2          ! Num_Weyls\n"
                "Direct     ! Direct or Cartesian coordinate\n"
                "0.0001       ! Radius of the ball surround a Weyl point\n"
                "-0.50000000    0.00000000    0.00000000       ! Positions of Weyl points, No. of lines should larger than Num_weyls\n"
                "-0.00000000    0.50000000   -0.00000000\n\n")

        if isloto == "T" and dt is not None and bc is not None:
            f.write("LOTO_DT\n"
                    f"{dt[0][0]:20.10f}    {dt[0][1]:20.10f}    {dt[0][2]:20.10f}\n"
                    f"{dt[1][0]:20.10f}    {dt[1][1]:20.10f}    {dt[1][2]:20.10f}\n"
                    f"{dt[2][0]:20.10f}    {dt[2][1]:20.10f}    {dt[2][2]:20.10f}\n\n")

            f.write("LOTO_BC\n")
            for i in range(num_atom):
                f.write(f"{bc[i][0][0]:20.10f}    {bc[i][0][1]:20.10f}    {bc[i][0][2]:20.10f}\n"
                        f"{bc[i][1][0]:20.10f}    {bc[i][1][1]:20.10f}    {bc[i][1][2]:20.10f}\n"
                        f"{bc[i][2][0]:20.10f}    {bc[i][2][1]:20.10f}    {bc[i][2][2]:20.10f}\n")
            f.write("LOTO_INTERPOLATE\n"
                    "11 11 11\n\n")

    print("Finished successfully!")


if __name__ == '__main__':
    main()
