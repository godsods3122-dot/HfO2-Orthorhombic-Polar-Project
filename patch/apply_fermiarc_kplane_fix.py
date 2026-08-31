#!/usr/bin/env python3
"""Simphony `fermiarc.f90` — KPLANE_SLAB 카드를 무시하고 full BZ 를 계산하는 문제.

## 증상

`SlabArc_calc = T` 로 등주파수 표면 스펙트럼을 뽑을 때, KPLANE_SLAB 에

```
KPLANE_SLAB
-0.30 -0.20      ! 시작점
 0.60  0.00      ! 첫 벡터
 0.00  0.40      ! 둘째 벡터
```

라고 좁은 창을 지정해도, 실제로 계산되는 영역은
`k_a ∈ [-0.30, 0.70]`, `k_b ∈ [-0.20, 0.80]` — 즉 시작점에서 시작하는
**full BZ 한 칸**이다.  `arc.dat_l` 의 kx, ky 폭을 재 보면 정확히 |b1|, |b2| 다.

출력에 경고가 찍히긴 한다:

```
WARNING : Your setting of KPLANE_SLAB has been modified because QPI
          calculation requires the information of the full BZ.
The first modified vector in QPI:   1.0000  0.0000
The second modified vector in QPI:  0.0000  1.0000
```

문제는 **QPI 를 켜지 않아도 항상** 이 변환이 일어난다는 것.  그래서

- 관심 영역만 조밀하게 보려던 계획이 조용히 어긋나고,
- 그림 축을 지정한 범위로 라벨링하면 좌표가 통째로 틀린다.
  (BZ 를 넘어간 부분이 되접히므로 대칭이 깨진 것처럼 보인다.)

## 원인

`src/fermiarc.f90` 의 ceiling/floor 블록이 조건 없이 실행된다:

```fortran
!> ceiling if K2D_vec are positive, floor if K2D_vec are negative
do i  = 1, 2
    if (K2D_vec1(i)>0) then
          K2D_vec_a(i)= ceiling(K2D_vec1(i))
    else
          K2D_vec_a(i)= floor(K2D_vec1(i))
    endif
    ...
```

`ceiling(0.6) = 1`, `ceiling(0.4) = 1` 이라 어떤 좁은 창을 줘도 1 이 된다.

## 고침

full BZ 가 실제로 필요한 것은 자기상관(joint DOS)을 쓰는 QPI 뿐이다.
`SlabQPI_kplane_calc` 가 켜졌을 때만 올림하고, 아니면 카드 값을 그대로 쓴다.
경고 문구도 실제로 변환이 일어날 때만 찍는다.

되돌리려면 `git checkout src/fermiarc.f90`.

적용:
```bash
cd <simphony 루트>          # src/ 가 보이는 곳
python3 apply_fermiarc_kplane_fix.py
cd src && make
```
`module.f90` 을 안 건드리므로 `make clean` 은 필요 없다.
"""
import os
import sys

OLD = """     !> ceiling if K2D_vec are positive, floor if K2D_vec are negative
     do i  = 1, 2
         if (K2D_vec1(i)>0) then
               K2D_vec_a(i)= ceiling(K2D_vec1(i))
         else
               K2D_vec_a(i)= floor(K2D_vec1(i))
         endif

         if (K2D_vec2(i)>0) then
            K2D_vec_b(i)= ceiling(K2D_vec2(i))
         else
            K2D_vec_b(i)= floor(K2D_vec2(i))
         endif
      enddo 
      if (cpuid==0) then
         write(stdout, '(a)')'WARNING : Your setting of KPLANE_SLAB has been modified because QPI calculation requires the information of the full BZ. '
         write(stdout, '((a, 2f8.4))')'The first modified vector in QPI: ', K2D_vec_a
         write(stdout, '((a, 2f8.4))')'The second modified vector in QPI: ', K2D_vec_b
      endif
"""

NEW = """     !> QPI (joint DOS) needs the full BZ, so there the KPLANE_SLAB vectors are
     !> rounded outwards to whole reciprocal lattice vectors.  For a plain
     !> SlabArc/SlabSpintexture run this rounding silently replaces the window
     !> the user asked for by a full BZ cell, so do it only when QPI is on.
     if (SlabQPI_kplane_calc) then
        do i  = 1, 2
            if (K2D_vec1(i)>0) then
                  K2D_vec_a(i)= ceiling(K2D_vec1(i))
            else
                  K2D_vec_a(i)= floor(K2D_vec1(i))
            endif

            if (K2D_vec2(i)>0) then
               K2D_vec_b(i)= ceiling(K2D_vec2(i))
            else
               K2D_vec_b(i)= floor(K2D_vec2(i))
            endif
         enddo
         if (cpuid==0) then
            write(stdout, '(a)')'WARNING : Your setting of KPLANE_SLAB has been modified because QPI calculation requires the information of the full BZ. '
            write(stdout, '((a, 2f8.4))')'The first modified vector in QPI: ', K2D_vec_a
            write(stdout, '((a, 2f8.4))')'The second modified vector in QPI: ', K2D_vec_b
         endif
      else
         K2D_vec_a= K2D_vec1
         K2D_vec_b= K2D_vec2
      endif
"""


def main():
    path = os.path.join('src', 'fermiarc.f90')
    if not os.path.isfile(path):
        sys.exit('src/fermiarc.f90 을 찾을 수 없다. simphony 루트에서 실행할 것.')
    s = open(path).read()
    if 'if (SlabQPI_kplane_calc) then' in s and 'K2D_vec_a= K2D_vec1' in s:
        print('이미 적용되어 있다.  변경 없음.')
        return
    if OLD not in s:
        sys.exit('원본 블록을 찾지 못했다. fermiarc.f90 이 이미 수정되었거나 버전이 다르다.')
    open(path, 'w').write(s.replace(OLD, NEW, 1))
    print('src/fermiarc.f90 수정 완료.')
    print('이제  cd src && make  로 다시 빌드할 것 (make clean 불필요).')


if __name__ == '__main__':
    main()
