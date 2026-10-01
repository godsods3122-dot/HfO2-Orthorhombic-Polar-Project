#!/usr/bin/env python3
"""Simphony `fermiarc.f90` — SlabArc 가 LO-TO 보정을 건너뛴다.

## 증상

`LOTO_correction = T`, `LOTO_method = 'phonopy'` 로 설정하고 PN.out 에도

```
 LOTO_correction                   :  T
 We found LOTO_DT card for LOTO correction
 We found LOTO_BC card for LOTO correction
```

가 찍히는데도, `SlabArc_calc = T` 로 얻은 `arc.dat_l/_r/_bulk` 는
LO-TO 가 **전혀 반영되지 않은** 표면 스펙트럼이다.

## 원인

같은 표면 Green 함수를 쓰는 두 루틴이 서로 다르게 짜여 있다.

`src/surfstat.f90` (SlabSS) 는 제대로 갈래를 탄다:

```fortran
if (index(Particle,'phonon')/=0.and.LOTO_correction) then
   call ham_qlayer2qlayer_LOTO(k,H00,H01)
else
   call ham_qlayer2qlayer(k,H00,H01)
endif
```

`src/fermiarc.f90` (SlabArc / SlabSpintexture / SlabQPI) 는 조건 없이:

```fortran
call ensure_LR_realspace()          ! 장거리항을 실공간에 만들어는 둔다
call ham_qlayer2qlayer(k,H00,H01)   ! ...그런데 안 쓰는 판을 부른다
```

`ham_qlayer2qlayer` 는 `HmnR_LRslab` 를 더하지 않고 per-k 비해석항도 넣지
않는다.  바로 위의 `ensure_LR_realspace()` 는 결과에 아무 영향이 없는
헛일이 된다.

`gapshape3D` 가 `ham_bulk_latticegauge` 를 무조건 부르던 것과 같은 종류의
누락이다 (`apply_gapcube_loto_fix.py` 참조).

## 고침

`surfstat.f90` 과 똑같은 갈래를 넣는다.

되돌리려면 `git checkout src/fermiarc.f90`.

적용:
```bash
cd <simphony 루트>          # src/ 가 보이는 곳
python3 apply_fermiarc_loto_fix.py
cd src && make
```
`module.f90` 을 안 건드리므로 `make clean` 은 필요 없다.
"""
import os
import sys

OLD = """        call ham_qlayer2qlayer(k,H00,H01)
        call now(time2)"""

NEW = """        !> deal with phonon system: the LO-TO corrected layer Hamiltonian is a
        !> different routine, exactly as surfstat.f90 branches.  Without this the
        !> ensure_LR_realspace() call above has no effect and the arc spectrum is
        !> computed with no LO-TO at all.
        if (index(Particle,'phonon')/=0.and.LOTO_correction) then
           call ham_qlayer2qlayer_LOTO(k,H00,H01)
        else
           call ham_qlayer2qlayer(k,H00,H01)
        endif
        call now(time2)"""


def main():
    path = os.path.join('src', 'fermiarc.f90')
    if not os.path.isfile(path):
        sys.exit('src/fermiarc.f90 을 찾을 수 없다. simphony 루트에서 실행할 것.')
    s = open(path).read()
    if 'call ham_qlayer2qlayer_LOTO(k,H00,H01)' in s:
        print('이미 적용되어 있다.  변경 없음.')
        return
    if s.count(OLD) != 1:
        sys.exit('원본 호출부를 유일하게 찾지 못했다 (%d 건). 버전이 다르다.'
                 % s.count(OLD))
    open(path, 'w').write(s.replace(OLD, NEW, 1))
    print('src/fermiarc.f90 수정 완료.')
    print('이제  cd src && make  로 다시 빌드할 것 (make clean 불필요).')


if __name__ == '__main__':
    main()
