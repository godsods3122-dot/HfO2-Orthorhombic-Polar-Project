#!/usr/bin/env python3
"""
phonopy2TBDAT.py

phonopy의 FORCE_SETS + POSCAR (조화 2차 힘 상수)를 읽어서
wannier90/wanniertools/phonnier가 요구하는 실공간 phonon
tight-binding 해밀토니안 파일(name_phononTB_hr.dat, wannier90_hr.dat
형식)을 생성한다.

이 스크립트는 QE 기준으로 작성된 QE2TBDAT.py의 phonopy 버전이다.
아래 두 가지를 QE 버전과 "물리적으로 1:1 대응"되도록 유지했다:

  1) smallest_vectors(): 원자쌍 사이의 최단(Wigner-Seitz) 이미지를
     찾고, 동일 거리에 여러 개의 등가 이미지(degenerate image)가
     있으면 그 개수(multi)로 나눠서 기여도를 균등하게 분배한다.
     이것이 바로 phonopy가 내부적으로 (NAC 없는) 동역학 행렬을
     Fourier 보간할 때 쓰는 것과 같은 알고리즘이며, wannier90의
     "use_ws_distance"와 동일한 사상이다. 이 분배를 제대로 하지
     않으면 고대칭점에서 축퇴가 깨지는 원인이 된다.

  2) [중요, 뒤늦게 수정됨] --born을 주면 phonopy 자신의
     DynamicalMatrixGL.make_Gonze_nac_dataset()으로 장거리
     dipole-dipole 성분이 빠진 단거리(short-range) 힘상수를 뽑아서
     사용한다. finite-displacement supercell 힘상수(ph.force_constants)
     는 "commensurate q점에서는 이미 dipole-dipole 상호작용이 포함된
     total 힘상수"이다(phonopy/phono3py 구현 논문, Euphonic 문서 참조) -
     극성 물질에서 이 상태 그대로 hr.dat을 만들면, Simphony가 LOTO_BC로
     analytic dipole-dipole 항을 다시 더할 때 이중 계산(double counting)이
     된다. 반드시 hr.dat 생성 "전" 단계에서 빼야 한다.

  2) write 단계에서 최종 R점들의 degeneracy는 전부 1로 쓴다
     (QE 버전의 write_phonon_hr/write_phonon_hr2와 동일한 관례).
     WS 축퇴는 이미 1)에서 multi로 흡수했기 때문이다.

QE 버전에서 실제로 "죽어 있던"(주석 처리된) write_phonon_hr 경로에는
분수좌표 벡터와 데카르트 좌표를 단위 구분 없이 더하는 버그가 있어
그대로 이식하지 않았다. 대신 좌표 변환을 전 구간에서 명시적으로
데카르트<->분수 좌표로 구분해 다시 구현했다. (질문하신 "축퇴가 깨지는"
현상의 유력한 원인 중 하나가 바로 이런 좌표계 혼동이다.)

사용법:
    python phonopy2TBDAT.py POSCAR FORCE_SETS "2 2 2" -o name_phononTB_hr.dat
    python phonopy2TBDAT.py POSCAR FORCE_SETS "2 0 0 0 2 0 0 0 2" -p "0.5 0.5 0 0.5 0 0.5 0 0.5 0.5"

supercell_matrix / --primitive-matrix 는
  - 공백으로 구분된 숫자 3개  -> 대각행렬
  - 공백으로 구분된 숫자 9개  -> 3x3 행렬 (row-major)
  - --primitive-matrix 는 "auto"/"identity" 같은 phonopy 키워드도 허용
"""

import argparse
from datetime import datetime

import numpy as np

from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_SETS, parse_BORN
from phonopy.interface.calculator import read_crystal_structure
from phonopy.structure.cells import get_reduced_bases

# phonopy의 잘 알려진 "VASP 힘상수 단위(eV/Angstrom^2) -> THz" 변환 계수.
# sqrt(eigenvalue[eV/AMU/Angstrom^2]) * VaspToTHz = 주파수[THz].
# 최신 phonopy에서는 phonopy.units가 deprecated 되었으므로 값을 직접
# 상수로 박아 넣어 경고/버전 의존성을 피한다.
VaspToTHz = 15.633302300230191


# ----------------------------------------------------------------------
# 유틸리티
# ----------------------------------------------------------------------

def printProgressBar(iteration, total, prefix="", suffix="", decimals=1,
                      length=30, fill="█", printEnd="\r"):
    if total == 0:
        return
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print(f"\r{prefix} |{bar}| {percent}% {suffix}", end=printEnd)
    if iteration == total:
        print()


def covariant_coordinates(basis, vectors):
    """
    데카르트 벡터(들)를 basis(행벡터가 격자벡터)의 성분으로 변환한다.
    vectors = c @ basis 를 만족하는 c를 반환.
    vectors: (...,3) 또는 (3,)
    """
    inv_basis = np.linalg.inv(basis)
    return np.dot(vectors, inv_basis)


def parse_matrix_arg(s):
    nums = [float(x) for x in s.replace(",", " ").split()]
    if len(nums) == 3:
        return np.diag(nums)
    if len(nums) == 9:
        return np.array(nums).reshape(3, 3)
    raise argparse.ArgumentTypeError(
        "행렬 인자는 숫자 3개(대각) 또는 9개(3x3, row-major)여야 합니다: '{}'".format(s)
    )


_BRAVAIS_KEYWORDS = {"p": "P", "f": "F", "i": "I", "a": "A", "c": "C", "r": "R"}


def resolve_nac_factor(user_value=None):
    """
    NAC(비해석항) 단위 변환 계수를 결정한다.

    VASP 기준(힘상수가 eV/Angstrom^2)에서 이 값은 Hartree * Bohr
    = 14.399652 eV*Angstrom 이다. phonopy 버전에 따라 물리 상수를 얻는
    API가 다르므로 순서대로 시도하고, 모두 실패하면 하드코딩 값을 쓴다.
    """
    if user_value is not None:
        return user_value
    try:  # phonopy v4+
        from phonopy.physical_units import get_physical_units
        u = get_physical_units()
        return u.Hartree * u.Bohr
    except Exception:
        pass
    try:  # phonopy v3 이하
        from phonopy.units import Hartree, Bohr
        return Hartree * Bohr
    except Exception:
        pass
    return 14.399652


def parse_primitive_matrix_arg(s):
    key = s.strip().lower()
    if key in ("identity", "none", "unit", "unitcell"):
        # POSCAR 자체가 이미 (FORCE_SETS 계산에 쓰인) 원시셀인 경우.
        # 주의: phonopy v4부터 primitive_matrix=None은 새 기본값 'auto'로
        # 풀린다(PrimitiveMatrixAutoDefaultWarning, v3까지는 None=identity
        # 였음). None이 아니라 명시적으로 'P'(Bravais: 센터링 없음=identity)
        # 를 넘겨서 v3/v4 어느 쪽이든 원시셀이 unit cell과 좌표계(방향)까지
        # 완전히 동일하게 유지되도록 강제한다.
        return "P"
    if key == "auto":
        # 주의: spglib 기반 자동 축소/표준화를 쓰며, 대칭이 낮거나(예: 변형된
        # biaxial-strain 구조) 원자좌표에 잡음이 있으면 spglib이 unit cell과
        # "다른 방향(회전/전단)"의 primitive 기저를 돌려줄 수 있다. 그 경우
        # supercell의 데카르트 벡터를 그 primitive 기저로 사영해도 정수가
        # 나오지 않아 이 스크립트가 에러를 낸다. 원래 POSCAR가 이미 primitive
        # 셀이라면 'identity'를 쓰는 것이 안전하다.
        return "auto"
    if key in _BRAVAIS_KEYWORDS:
        return _BRAVAIS_KEYWORDS[key]
    return parse_matrix_arg(s)


# ----------------------------------------------------------------------
# 기하학: 최단 이미지 벡터 탐색 (Wigner-Seitz, 축퇴 포함)
# ----------------------------------------------------------------------

def get_smallest_vectors(pos_from_frac, pos_to_frac, sbasis, symprec=1e-5):
    """
    pos_from_frac: (n_from, 3) supercell 기준 분수좌표
    pos_to_frac  : (n_to, 3)   supercell 기준 분수좌표
    sbasis       : (3,3) supercell 격자벡터 (행벡터, 데카르트)

    각 (i in n_from, j in n_to) 쌍에 대해 i -> j 로 향하는 최단
    이미지 벡터(들)를 데카르트 좌표로 반환한다. 동일 최단거리를 갖는
    이미지가 여러 개면 전부 반환하고 개수를 multi에 기록한다.

    반환
    ----
    svecs : (n_from, n_to, 27, 3) 데카르트 벡터 (0-padding)
    multi : (n_from, n_to) int, 유효한(축퇴) 이미지 개수
    """
    n_from = len(pos_from_frac)
    n_to = len(pos_to_frac)

    # phonopy 방식: (0,a,b,c,-a-b-c,-a,-b,-c,a+b+c) 조합으로 만든
    # 65개(중복 제거 후)의 격자점 후보 안에서 최근접 이미지를 찾는다.
    lattice_1D = (-1, 0, 1)
    lattice_4D = np.array(
        [[i, j, k, l] for i in lattice_1D for j in lattice_1D
         for k in lattice_1D for l in lattice_1D]
    )
    bases = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, -1, -1]])
    lattice_points = np.unique(np.dot(lattice_4D, bases), axis=0)
    npoints = len(lattice_points)

    svecs = np.zeros((n_from, n_to, 27, 3))
    multi = np.zeros((n_from, n_to), dtype=np.int32)

    total = n_from * n_to
    counter = 0
    printProgressBar(counter, total, prefix="최단벡터 탐색:")
    for i in range(n_from):
        for j in range(n_to):
            # i -> j 로 향하는 벡터 (분수좌표계에서 격자점 후보를 더함)
            frac_diff = pos_to_frac[j] - pos_from_frac[i] + lattice_points
            cart_diff = frac_diff @ sbasis
            length = np.linalg.norm(cart_diff, axis=1)
            minimum = length.min()
            close = np.where(np.abs(length - minimum) < symprec)[0]
            if len(close) > 27:
                raise RuntimeError(
                    "원자쌍 ({}, {})에서 27개보다 많은 등거리 이미지가 발견되었습니다. "
                    "구조/symprec을 확인하세요.".format(i, j)
                )
            svecs[i, j, : len(close)] = cart_diff[close]
            multi[i, j] = len(close)
            counter += 1
            printProgressBar(counter, total, prefix="최단벡터 탐색:")

    return svecs, multi


# ----------------------------------------------------------------------
# 실공간 해밀토니안 조립 + wannier90_hr.dat 형식 출력
# ----------------------------------------------------------------------

def build_and_write_hr(fc, svecs, multi, masses, p2s_map, s2p_map, p2p_map,
                        fc_row_index,
                        pos_cart, num_patom, num_satom, primitive_cell, outfile,
                        factor, sdim=20, tol_int=1e-3, min_amp=1e-20):
    """
    fc : (num_satom, num_satom, 3, 3) phonopy full force constants
    svecs, multi : get_smallest_vectors()의 결과, shape (num_patom, num_satom, ...)
                   (pos_from = 원시셀 원자들의 supercell 내 홈 이미지)
    masses : (num_patom,) 원시셀 질량
    p2s_map, s2p_map, p2p_map : phonopy Primitive의 매핑
    pos_cart : (num_satom, 3) supercell 원자들의 데카르트 좌표
    primitive_cell : (3,3) 원시셀 격자벡터 (행벡터, 데카르트)
    factor : 주파수 단위 변환 계수 (예: phonopy.units.VaspToTHz).
             H(R) 블록에 factor**2를 곱해서, 나중에 대각화한 뒤
             sqrt(eigenvalue)가 곧바로 THz가 되도록 한다.
    sdim : 원시격자 R벡터 성분의 최대 절대값 한계 (이 범위를 넘으면 에러)

    핵심 포인트
    -----------
    svecs[iatom, katom, rr]는 "원자 i(홈 이미지) -> 원자 katom" 최단
    분리벡터(raw separation)이며, 여기에는 순수 격자 변환 R뿐 아니라
    두 원자의 unit-cell 내부 부격자(sublattice) 위치 차이
    (d_katom_home - d_i_home)까지 섞여 있다. 이 부격자 위치차는 일반적으로
    원시격자 벡터의 정수배가 아니므로, raw separation을 곧바로 반올림하면
    안 된다. home_i, home_j(각각 p2s_map[iatom], s2p_map[katom]의 데카르트
    위치)를 더하고 빼서 순수 격자벡터 R만 분리한 뒤에 반올림해야 한다.
    """
    norbs = num_patom * 3
    dim = 2 * sdim + 1
    hr_mat = np.zeros((dim, dim, dim, norbs, norbs), dtype=np.complex128)

    total = num_patom * num_satom
    counter = 0
    printProgressBar(counter, total, prefix="H(R) 조립:")
    for iatom in range(num_patom):
        i_s = p2s_map[iatom]        # 기하학(위치)용: supercell 인덱스
        fc_row = fc_row_index[iatom]  # 힘상수 배열 인덱싱용 (compact/full 구분)
        home_i = pos_cart[i_s]
        for katom in range(num_satom):
            rep_katom = s2p_map[katom]      # katom의 부격자를 대표하는 홈 이미지의 supercell 인덱스
            jatom = p2p_map[rep_katom]      # 원시셀 원자 인덱스 (0..num_patom-1)
            home_j = pos_cart[rep_katom]
            m = int(multi[iatom, katom])
            if m == 0:
                counter += 1
                continue
            mass_sqrt = np.sqrt(masses[iatom] * masses[jatom])
            for rr in range(m):
                R_raw = svecs[iatom, katom, rr]           # i(홈) -> katom 최단 raw 분리벡터
                R_cart = R_raw + home_i - home_j           # 부격자 위치차 제거 -> 순수 격자벡터
                R_frac = covariant_coordinates(primitive_cell, R_cart)
                R_int = np.round(R_frac).astype(int)
                if np.max(np.abs(R_frac - R_int)) > tol_int:
                    raise RuntimeError(
                        "원자 {} <-> supercell 원자 {} 사이의 격자벡터가 원시격자의 "
                        "정수 배가 아닙니다 (R_frac={}). supercell이 원시셀과 "
                        "commensurate 하지 않거나 symprec/tol_int이 부적절할 수 "
                        "있습니다.".format(iatom, katom, R_frac)
                    )
                rx, ry, rz = R_int + sdim
                if not (0 <= rx < dim and 0 <= ry < dim and 0 <= rz < dim):
                    raise RuntimeError(
                        "R벡터가 허용 범위(sdim={})를 벗어났습니다: {}. "
                        "--sdim 값을 늘려서 다시 시도하세요.".format(sdim, R_int)
                    )
                block = (fc[fc_row, katom] / (m * mass_sqrt)) * (factor ** 2)
                hr_mat[rx, ry, rz,
                       iatom * 3: iatom * 3 + 3,
                       jatom * 3: jatom * 3 + 3] += block
            counter += 1
            printProgressBar(counter, total, prefix="H(R) 조립:")

    # 0이 아닌 R점만 모으기
    rpts = []
    blocks = []
    for rx in range(dim):
        for ry in range(dim):
            for rz in range(dim):
                blk = hr_mat[rx, ry, rz]
                if np.abs(blk).sum() < min_amp:
                    continue
                rpts.append((rx - sdim, ry - sdim, rz - sdim))
                blocks.append(blk)
    nrpt = len(rpts)
    dege_rpts = np.ones(nrpt, dtype=np.int32)

    with open(outfile, "w") as f:
        f.write(" Written by phonopy2TBDAT.py on " + str(datetime.now()) + "\n")
        f.write("{:>12d}\n".format(norbs))
        f.write("{:>12d}\n".format(nrpt))
        nl = int(np.ceil(nrpt / 15.0)) if nrpt > 0 else 0
        for n in range(nl):
            line = "    " + "    ".join(
                str(int(d)) for d in dege_rpts[n * 15:(n + 1) * 15]
            )
            f.write(line + "\n")
        for irpt in range(nrpt):
            rx, ry, rz = rpts[irpt]
            block = blocks[irpt]
            for jorb in range(norbs):
                for iorb in range(norbs):
                    rp = block[iorb, jorb].real
                    ip = block[iorb, jorb].imag
                    f.write(
                        "{:8d}{:8d}{:8d}{:8d}{:8d}{:20.10f}{:20.10f}\n".format(
                            rx, ry, rz, iorb + 1, jorb + 1, rp, ip
                        )
                    )

    return nrpt, norbs


# ----------------------------------------------------------------------
# 메인
# ----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="phonopy FORCE_SETS + POSCAR -> wannier90_hr.dat 형식 phonon TB 해밀토니안 생성"
    )
    parser.add_argument("poscar", help="phonopy 원시(또는 unit) 셀 POSCAR 경로")
    parser.add_argument("force_sets", help="FORCE_SETS 경로")
    parser.add_argument(
        "supercell_matrix", type=parse_matrix_arg,
        help="supercell 행렬. '2 2 2'(대각) 또는 9개 숫자(3x3, row-major)"
    )
    parser.add_argument(
        "--born", type=str, default=None,
        help="BORN 파일 경로 (극성 물질, LOTO_correction을 쓸 경우 필수). "
             "주면 phonopy의 Gonze-Lee NAC 기능으로 dipole-dipole 성분이 "
             "제거된 단거리 힘상수를 hr.dat에 쓴다. input.py에도 반드시 "
             "동일한 --born/-p 조합을 넘겨서 LOTO_DT/LOTO_BC와 짝을 맞출 것."
    )
    parser.add_argument(
        "-p", "--primitive-matrix", default="identity", type=parse_primitive_matrix_arg,
        help="primitive matrix. 'identity'(POSCAR 자체가 원시셀, 기본값) / "
             "'auto'(spglib으로 자동 축소, 좌표계가 회전될 수 있어 주의) / "
             "'P','F','I','A','C','R' / '3 또는 9개 숫자'"
    )
    parser.add_argument("-o", "--outfile", default="name_phononTB_hr.dat",
                         help="출력 파일 이름 (기본값: name_phononTB_hr.dat)")
    parser.add_argument("--symprec", type=float, default=1e-5,
                         help="최단벡터 탐색 시 축퇴 판정 허용오차 (기본값: 1e-5)")
    parser.add_argument("--sdim", type=int, default=20,
                         help="원시격자 R벡터 성분의 최대 절대값 한계 (기본값: 20)")
    parser.add_argument("--factor", type=float, default=None,
                         help="주파수 변환 계수(기본값: phonopy VaspToTHz). "
                              "FORCE_SETS가 VASP 단위(eV/Angstrom^2)가 아니라면 지정하세요.")
    parser.add_argument(
        "--nac-factor", type=float, default=None,
        help="NAC 단위 변환 계수. 지정하지 않으면 BORN 파일 값을 쓰고, "
             "거기에도 없으면 VASP 기준 Hartree*Bohr(=14.399652 eV*Angstrom)를 "
             "사용한다. FORCE_SETS가 eV/Angstrom^2 단위가 아니면 지정하세요."
    )
    parser.add_argument(
        "--fc-symmetrize-level", type=int, default=3,
        help="힘상수에 translational invariance(음향 총합 규칙, ASR) + "
             "permutation symmetry를 강제하는 반복 횟수 (기본값: 3). "
             "finite-displacement로 얻은 FC는 기본적으로 ASR을 정확히 "
             "만족하지 않아 Gamma 근처 acoustic 브랜치가 미세하게 허수로 "
             "나올 수 있다 - 이를 교정한다. 0을 주면 끌 수 있다."
    )
    parser.add_argument(
        "--fc-spg-symmetrize", dest="fc_spg_symmetrize", action="store_true", default=True,
        help="힘상수에 결정의 공간군(space group) 대칭을 강제한다 (기본값: 켜짐). "
             "DFT 잡음으로 인해 대칭으로 보호되는 밴드 교차/분기점에서 "
             "스퓨리어스하게 갈라지거나 휘는 문제를 교정한다."
    )
    parser.add_argument(
        "--no-fc-spg-symmetrize", dest="fc_spg_symmetrize", action="store_false",
        help="공간군 대칭 강제를 끈다."
    )
    args = parser.parse_args()

    factor = args.factor if args.factor is not None else VaspToTHz

    print("입력 파일을 읽는 중입니다...")
    unitcell, _ = read_crystal_structure(filename=args.poscar, interface_mode="vasp")

    ph = Phonopy(
        unitcell,
        supercell_matrix=args.supercell_matrix,
        primitive_matrix=args.primitive_matrix,
    )

    force_sets = parse_FORCE_SETS(filename=args.force_sets)
    ph.dataset = force_sets
    ph.produce_force_constants(calculate_full_force_constants=True)

    if args.fc_spg_symmetrize:
        # 결정의 공간군(space group) 대칭을 힘상수에 강제한다. finite-
        # displacement로 얻은 FC는 DFT 잡음 때문에 결정 자체의 대칭을
        # 완벽히 만족하지 못하는 경우가 흔한데, 이 잡음은 대칭으로 보호되는
        # 밴드 교차/분기점에서만 스퓨리어스하게 갈라지거나 휘는 형태로
        # 나타난다 (그 외 구간에서는 잘 안 보임). WT 공식 phonon_hr.py도
        # 이 단계를 ASR/permutation 대칭화보다 먼저 적용한다.
        print("결정의 공간군 대칭을 힘상수에 강제합니다 (spglib 기반)...")
        ph.symmetrize_force_constants_by_space_group()

    # ASR(음향 총합 규칙) 위반 정도 진단: 각 원자 i에 대해 sum_j FC[i,j]가
    # 0이어야 하는데, finite-displacement로 얻은 FC는 보통 정확히 0이 아니다.
    drift_before = np.max(np.abs(np.sum(ph.force_constants, axis=1)))
    print("교정 전 ASR 위반 정도(최대 drift): {:.6e}".format(drift_before))

    if args.fc_symmetrize_level > 0:
        print("힘상수에 음향 총합 규칙(ASR)/permutation symmetry를 "
              "강제합니다 (level={})...".format(args.fc_symmetrize_level))
        ph.symmetrize_force_constants(level=args.fc_symmetrize_level)
        drift_after = np.max(np.abs(np.sum(ph.force_constants, axis=1)))
        print("교정 후 ASR 위반 정도(최대 drift): {:.6e}".format(drift_after))

    if args.born:
        # finite-displacement supercell 힘상수는 commensurate q점에서는
        # 이미 dipole-dipole(장거리) 상호작용을 포함한 "total" 힘상수다.
        # Simphony(LOTO_correction=T)는 hr.dat이 이게 빠진 short-range만
        # 담고 있다고 가정하고 LOTO_BC로 analytic dipole-dipole 항을
        # 별도로 더하므로, 여기서 미리 빼지 않으면 이중 계산된다.
        # phonopy 자신의 Gonze-Lee NAC 구현을 그대로 재사용해 이 subtraction을
        # 수행한다 (직접 Ewald 합을 재구현하는 것보다 안전).
        print("BORN 파일을 읽어 단거리(dipole-dipole 제거) 힘상수를 "
              "추출합니다 (phonopy Gonze-Lee NAC)...")
        nac_params = parse_BORN(ph.primitive, filename=args.born)
        # NAC 단위 변환 계수. BORN 파일 첫 줄에 값이 적혀 있으면 parse_BORN이
        # 채워주지만, 없으면 'factor' 키 자체가 없어서 KeyError가 난다.
        # (phonopy CLI는 계산기 종류에서 기본값을 자동으로 채우지만 API는 안 채움)
        # VASP 기준 값은 Hartree * Bohr = 14.399652 eV*Angstrom.
        if "factor" not in nac_params:
            nac_params["factor"] = resolve_nac_factor(args.nac_factor)
            print("  BORN 파일에 NAC 계수가 없어 기본값을 사용합니다: "
                  "factor = {:.6f}".format(nac_params["factor"]))
        else:
            print("  BORN 파일에 적힌 NAC 계수를 사용합니다: "
                  "factor = {:.6f}".format(nac_params["factor"]))
        # NAC method를 Gonze로 명시한다. (phonopy 기본값도 GONZE지만,
        # BORN 파일이나 버전에 따라 달라질 수 있으므로 강제한다.)
        nac_params["method"] = "gonze"
        ph.nac_params = nac_params

        total_fc_before = ph.force_constants.copy()

        # Gonze-Lee DynamicalMatrix(DynamicalMatrixGL)를 실제로 생성시키기
        # 위해 밴드 계산을 1회 트리거한다. 결과 밴드 자체는 쓰지 않고
        # dm 객체만 필요하다.
        ph.run_band_structure([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]],
                               is_band_connection=False)
        dm = ph._dynamical_matrix
        if not hasattr(dm, "make_Gonze_nac_dataset"):
            raise RuntimeError(
                "동역학 행렬이 Gonze-Lee(DynamicalMatrixGL) 타입이 아닙니다 "
                "(실제 타입: {}). BORN 파일과 phonopy 버전을 확인하세요."
                .format(type(dm).__name__)
            )
        dm.make_Gonze_nac_dataset()
        srfc = dm.short_range_force_constants

        dd_removed = np.max(np.abs(total_fc_before - srfc))
        print("제거된 dipole-dipole 성분의 최대 크기: {:.6e}".format(dd_removed))
        drift_sr = np.max(np.abs(np.sum(srfc, axis=1)))
        print("단거리 추출 후 ASR 위반 정도(최대 drift): {:.6e}".format(drift_sr))

        # 중요: SR 힘상수를 원래 ph(=nac_params가 붙어있는 객체)에 덮어쓰면
        # 안 된다. force_constants를 대입하는 순간 캐시된 동역학 행렬이
        # 무효화되고, nac_params가 살아있으면 이후 NAC이 다시 얹힌 상태로
        # 재구성될 수 있다(=dipole-dipole 이중 계산). 레퍼런스 구현이
        # load_phonopy(None)로 별도 객체를 만드는 이유가 이것이다.
        # 여기서도 NAC이 전혀 없는 깨끗한 Phonopy 객체를 새로 만들어
        # 거기에 SR 힘상수만 넣고, 이후 hr.dat 조립에 그 객체를 쓴다.
        ph_sr = Phonopy(
            unitcell,
            supercell_matrix=args.supercell_matrix,
            primitive_matrix=args.primitive_matrix,
        )
        ph_sr.force_constants = srfc
        ph = ph_sr
    else:
        print("*** 경고: --born 미지정 -> dipole-dipole 뺄셈이 전혀 "
              "적용되지 않습니다. hr.dat은 total 힘상수 그대로입니다. "
              "극성 물질에서 Simphony의 LOTO_correction=T와 함께 쓰면 "
              "이중 계산됩니다. ***")

    fc = ph.force_constants  # compact (num_patom, num_satom, 3, 3) 또는
                             # full (num_satom, num_satom, 3, 3)
    primitive = ph.primitive
    supercell = ph.supercell

    masses = primitive.masses
    num_patom = len(primitive)
    num_satom = len(supercell)

    p2s_map = primitive.p2s_map
    s2p_map = primitive.s2p_map
    p2p_map = primitive.p2p_map

    sbasis = supercell.cell        # (3,3) 데카르트, supercell 격자벡터 (원본)
    pcell = primitive.cell         # (3,3) 데카르트, 원시셀 격자벡터

    spos_all = supercell.scaled_positions
    pos_cart = spos_all @ sbasis   # supercell 원자들의 데카르트 좌표 (기저 선택과 무관, 항상 정확)

    # WT(WannierTools) 공식 phonon_hr.py와 동일하게, Wigner-Seitz 최단벡터
    # 탐색은 supercell의 원래 기저가 아니라 격자 축소(Niggli reduction)된
    # 기저에서 수행한다. phonopy 자신도 내부적으로 이 축소된 기저를 써서
    # 최단벡터/축퇴도를 계산하는데, 대칭성이 낮은 구조(예: orthorhombic)에서는
    # 원래 기저 그대로 탐색할 때와 어떤 이미지들이 "동일 거리"로 축퇴되는지
    # 판정이 달라질 수 있다. 여기서 그 차이를 없앤다.
    reduced_bases = get_reduced_bases(sbasis, tolerance=args.symprec)
    pos_frac_reduced = pos_cart @ np.linalg.inv(reduced_bases)
    pos_from_reduced = pos_frac_reduced[p2s_map]

    print("입력 완료. 원시셀 원자수={}, supercell 원자수={}".format(num_patom, num_satom))

    # 힘상수 배열 형식 판별. phonopy는 두 가지 형식을 쓴다:
    #   full   : (num_satom, num_satom, 3, 3)  -> 첫 축을 supercell 인덱스(p2s_map)로 접근
    #   compact: (num_patom, num_satom, 3, 3)  -> 첫 축을 원시셀 인덱스(0..num_patom-1)로 접근
    # Gonze 단거리 힘상수(short_range_force_constants)는 compact로 나올 수 있는데,
    # 이때 full 형식처럼 p2s_map으로 인덱싱하면 IndexError 없이 "조용히" 엉뚱한
    # 행을 읽어가므로(예: p2s_map=[0,8,16,...]) 반드시 구분해야 한다.
    print("힘상수 배열 shape: {}".format(fc.shape))
    if fc.shape[0] == num_satom:
        print("  -> full 형식으로 인식 (첫 축 = supercell 인덱스)")
        fc_row_index = np.asarray(p2s_map)
    elif fc.shape[0] == num_patom:
        print("  -> compact 형식으로 인식 (첫 축 = 원시셀 인덱스)")
        fc_row_index = np.arange(num_patom)
    else:
        raise RuntimeError(
            "힘상수 배열의 첫 축 크기({})가 supercell 원자수({})도 "
            "원시셀 원자수({})도 아닙니다. 형식을 판별할 수 없습니다."
            .format(fc.shape[0], num_satom, num_patom)
        )

    print("Wigner-Seitz 최단벡터 및 축퇴도를 계산합니다 (격자 축소 기저 사용)...")
    svecs, multi = get_smallest_vectors(
        pos_from_reduced, pos_frac_reduced, reduced_bases, symprec=args.symprec
    )

    print("실공간 해밀토니안을 조립하고 {} 파일을 씁니다...".format(args.outfile))
    nrpt, norbs = build_and_write_hr(
        fc, svecs, multi, masses, p2s_map, s2p_map, p2p_map,
        fc_row_index,
        pos_cart, num_patom, num_satom, pcell, args.outfile,
        factor=factor, sdim=args.sdim,
    )

    print("완료. norbs = {}, nrpt = {}".format(norbs, nrpt))


if __name__ == "__main__":
    main()
