#!/usr/bin/env python3
"""
apply_gapcube_loto_fix.py — BulkGap_cube_calc 가 LO-TO 보정을 건너뛰는 버그를 고친다.

기존 `apply_simphony_loto_fix.py` 의 후속 패치다. 그 패치는 `ek_bulk`, `ek_slab`,
`wanniercenter` 의 LO-TO 경로를 고쳤지만, `fermisurface.f90` 의 `gapshape3D`
(= `BulkGap_cube_calc`) 는 손대지 않았다. 그 루틴은

    call ham_bulk_latticegauge(k, Hamk_bulk)

를 **무조건** 부른다. 극성 물질에서는 hr.dat 에 dipole-dipole 성분이 빠져 있으므로
(phonopy2TBDAT --born 이 빼서 쓴다), 이 경로는 **단거리 힘상수만의 스펙트럼**을 준다.
그 결과 BulkGap 이 내놓는 Ev/Ec/gap 이 실제 포논 밴드와 무관해지고, 노드 탐색이
조용히 엉뚱한 답을 낸다.

고침: `ek_bulk.f90` / `wanniercenter.f90` 과 동일하게 분기시킨다.
3D 격자에는 경로 방향이 없으므로 `LOTO_qdir_run` 은 설정하지 않는다 —
원 패치가 "3D grids, WCC 로 방향이 새지 않게 즉시 해제한다"고 명시한 용법 그대로다.

부수 사항 (고치지 않음, 알고 쓰면 됨): `gapshape3D` 가 쓰는 Ev/Ec 는 THz 주파수가
아니라 **동역학 행렬 고윳값** ω²(THz²) × 0.036749 (eV→Hartree) 다. 검증:
  Ev = 3.69304  <->  phonopy band17 = 10.02461 THz,  10.02461² × 0.036749 = 3.69304
따라서 파일의 gap 은 Δω 가 아니라 Δ(ω²) 이며, 노드 근방에서
  gap_file ≈ 0.073498 × ω̄(THz) × Δω(THz)
이다. 0 은 어느 쪽이든 0 이라 **노드 탐색에는 그대로 써도 된다.**

사용법:
    cd <simphony 루트>          # src/ 가 보이는 곳
    python3 apply_gapcube_loto_fix.py
    cd src && make clean && make
"""

import os
import sys

OLD = "         call ham_bulk_latticegauge(k, Hamk_bulk)\n"

NEW = """         !> deal with phonon system: gapshape3D used to bypass the LO-TO
         !> correction entirely, so BulkGap_cube_calc reported the gap of the
         !> dipole-subtracted short-range force constants. Branch the same way
         !> ek_bulk and wanniercenter already do. No LOTO_qdir_run is supplied:
         !> a 3D grid has no path direction, which is the documented usage.
         if (index(Particle,'phonon')/=0.and.LOTO_correction) then
            call ham_bulk_LOTO(k, Hamk_bulk)
         else
            call ham_bulk_latticegauge(k, Hamk_bulk)
         endif
"""

MARKER = "gapshape3D used to bypass the LO-TO"


def main():
    path = os.path.join("src", "fermisurface.f90")
    if not os.path.isfile(path):
        print("[!] src/fermisurface.f90 을 찾을 수 없습니다. simphony 루트에서 실행하세요.")
        sys.exit(1)

    text = open(path).read()
    if MARKER in text:
        print("[=] 이미 적용되어 있습니다. 아무것도 하지 않습니다.")
        return

    if "subroutine gapshape3D" not in text:
        print("[!] gapshape3D 를 찾을 수 없습니다. Simphony 버전을 확인하세요.")
        sys.exit(1)

    # gapshape3D 본문 안의 호출만 바꾼다 (같은 문장이 다른 루틴에도 있다).
    start = text.index("subroutine gapshape3D")
    end = text.index("end subroutine gapshape3D")
    body = text[start:end]
    if body.count(OLD) != 1:
        print("[!] gapshape3D 안에서 대상 호출을 %d 번 찾았습니다 (1 개를 기대). "
              "수동으로 확인하세요." % body.count(OLD))
        sys.exit(1)

    open(path, "w").write(text[:start] + body.replace(OLD, NEW) + text[end:])
    print("[+] src/fermisurface.f90 의 gapshape3D 를 패치했습니다.")
    print("    이어서: cd src && make clean && make")


if __name__ == "__main__":
    main()
