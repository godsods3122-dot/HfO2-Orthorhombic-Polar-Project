#!/usr/bin/env python3
"""Simphony `surfstat.f90` — `dos_l_only` / `dos_r_only` 가 항상 0 이 된다.

## 증상

`dos.dat_l` 의 4번째 열 `dos_l_only` (표면 상태만 남긴 값) 가 전 구간에서
`0.99999997E-09` (= eps9) 이다.  즉 "표면 전용" 출력이 쓸모가 없다.

## 원인

두 항의 궤도 개수가 다르다.

```fortran
do i= 1, NtopOrbitals            ! = Num_wann (한 단위포)
   io= TopOrbitals(i)
   dos_l(ikp, j)= dos_l(ikp,j)- aimag(GLL(io,io))
enddo
do i= 1, Ndim                    ! = Np * Num_wann (주층 전체)
   dos_bulk(ikp, j)= dos_bulk(ikp,j)- aimag(GB(i,i))
enddo
```

`dos_l` 은 표면 단위포 하나의 궤도 합, `dos_bulk` 은 주층(principal layer)
전체 = `Np` 개 단위포의 합이다.  그런데 뺄 때는 그냥 뺀다:

```fortran
dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
if (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
```

`Np = 2` 면 벌크를 두 배로 빼는 셈이라 결과가 **항상 음수**가 되고 전부
`eps9` 로 잘린다.  실제 HfO2 슬랩에서 확인: `dos_l - dos_bulk` 는 격자
100 %에서 음수, `dos_l - dos_bulk/Np` 로 바꾸면 30 %가 양수가 되면서
표면 상태 영역이 제대로 드러난다.

## 고침

궤도 개수 비로 맞춘다.  `NtopOrbitals = Num_wann` 인 보통의 경우
`dos_bulk/Np` 와 같다.

```fortran
dos_l_only = dos_l - dos_bulk * NtopOrbitals   / dble(Ndim)
dos_r_only = dos_r - dos_bulk * NBottomOrbitals/ dble(Ndim)
```

`fermiarc.f90` (SlabArc) 에는 `_only` 출력 자체가 없다.  같은 보정을
`arc.dat_l` 과 `arc.dat_bulk` 로 후처리에서 하면 된다
(`scripts/figs/fig3_arc.py` 의 `surface_only()` 참조).

되돌리려면 `git checkout src/surfstat.f90`.

적용:
```bash
cd <simphony 루트>
python3 apply_surfdos_only_norm_fix.py
cd src && make
```
"""
import os
import sys

OLD = """           dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
           dos_r_only(ikp, j)= dos_r_mpi(ikp, j)- dos_bulk_mpi(ikp, j)
           if (dos_r_only(ikp, j)<0) dos_r_only(ikp, j)=eps9"""

NEW = """           !> dos_l/dos_r sum over one surface unit cell (NtopOrbitals /
           !> NBottomOrbitals), while dos_bulk sums over the whole principal
           !> layer (Ndim = Np*Num_wann).  Rescale before subtracting, otherwise
           !> the bulk is over-subtracted by a factor Np and every _only value
           !> is clipped to eps9.
           dos_l_only(ikp, j)= dos_l_mpi(ikp, j)- dos_bulk_mpi(ikp, j) &
                               * dble(NtopOrbitals)/dble(Ndim)
           if (dos_l_only(ikp, j)<0) dos_l_only(ikp, j)=eps9
           dos_r_only(ikp, j)= dos_r_mpi(ikp, j)- dos_bulk_mpi(ikp, j) &
                               * dble(NBottomOrbitals)/dble(Ndim)
           if (dos_r_only(ikp, j)<0) dos_r_only(ikp, j)=eps9"""


def main():
    path = os.path.join('src', 'surfstat.f90')
    if not os.path.isfile(path):
        sys.exit('src/surfstat.f90 을 찾을 수 없다. simphony 루트에서 실행할 것.')
    s = open(path).read()
    if 'dble(NtopOrbitals)/dble(Ndim)' in s:
        print('이미 적용되어 있다.  변경 없음.')
        return
    if s.count(OLD) != 1:
        sys.exit('원본 블록을 유일하게 찾지 못했다 (%d 건). 버전이 다르다.' % s.count(OLD))
    open(path, 'w').write(s.replace(OLD, NEW, 1))
    print('src/surfstat.f90 수정 완료.')
    print('이제  cd src && make  로 다시 빌드할 것.')


if __name__ == '__main__':
    main()
