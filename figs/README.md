# PPT 그림 목록

모든 그림은 300 dpi PNG. 축 규약은 **표준 (a, b, c) — 편극축 = c**로 통일했다
(구조마다 native 축 순서가 달라서 `scripts/figs/axes.py` 가 경로를 매핑한다).
POSCAR 자체는 절대 치환하지 않는다 — FORCE_SETS 와 어긋난다.

| 파일 | 내용 | 만드는 법 |
|---|---|---|
| `fig1_dispersion_parent_pristine.png` | parent_pristine 포논 분산, Γ-X-S-Y-Γ-Z-U-R-T-Z | `run_band.py` → `fig1_2_bands.py` |
| `fig2_bands17_18_node.png` | 밴드 17/18 강조 + 노드 주변 절단 2장 | `fig2_node.py` |
| `fig3_fermi_arc.png` | ω = 10.0869 THz 등주파수 표면 스펙트럼 (bottom / top / bulk + 아크 분리) | `run_slab.py --mode arc` → `fig3_arc.py` |
| `fig4_slab_arc.png` | 두 사영점을 관통하는 직선 `k_a`=0.14649 위의 표면 스펙트럼 | `run_slab.py --mode ss` → `fig4_slabss.py` |
| `fig5_weyl_cone.png` | 바일 콘 3D + 주축 절단 (type-II 증거) | `fig5_cone.py` |
| `fig6_wcc_4weyl.png` | 4개 노드 각각의 구면 Wilson loop (WCC) | Simphony `WeylChirality_calc` → `fig6_8.py` |
| `fig7_strain_overlay.png` | biaxial strain 0 / −0.8 / −1.0 / −2.5 / −3.0 % 밴드 중첩 | `run_band.py` ×5 → `fig7_overlay.py` |
| `fig8_weyl_shift.png` | 4점 궤도의 이동 (unstrained → −0.8 %) | `fig6_8.py` |

## 숫자 요약 (발표에 그대로 쓸 수 있는 값)

| 항목 | 값 | 출처 |
|---|---|---|
| 바일 4점 위치 (k_c=0 평면) | `(±0.14649, ±0.07085, 0)` | Simphony BulkGap 401² + 정밀화 |
| 주파수 | 10.0869 THz = 41.72 meV | 〃 |
| chirality | `+1, −1, −1, +1` | Simphony `WeylChirality_calc` (PN.out) |
| gap (노드에서) | 6.5e-8 THz | phonopy |
| tilt 파라미터 T | 2.748 → **type-II** | 3D 국소 전개 |
| 면내 콘 기울기 | 1.1192 / 0.7619 THz/rlu | 〃 |
| 편극축(c) 방향 속도 | ≈0.0032 THz/rlu (면내의 1/350) | 〃 |
| −0.8 % 에서 위치 | `(±0.12942, ±0.06377, 0)` | 〃 |
| 이동량 | `Δk_a = −0.0171, Δk_b = −0.0071, |Δk| = 0.0185` | 〃 |
| 소멸 구간 | −0.8 % 와 −1.0 % 사이 | BulkGap + phonopy 둘 다 |

## 슬랩 설정 (fig 3, 4)

```
SURFACE  (1 0 0) / (0 1 0)     → 표면 법선 ∥ c = 편극축
NP = 2                          hr.dat 의 R3 범위가 ±1 이라 2 면 정확하다
LOTO_method = 'phonopy'         (빠뜨리면 밴드가 0.03 THz 조용히 어긋난다)
NumOccupied = 17,  SELECTED_OCCUPIED_BANDS 1-17
```

표면 법선을 c 로 잡아야 4개 사영점이 서로 겹치지 않는다. a 나 b 를 법선으로
잡으면 부호가 반대인 두 점이 같은 자리에 사영되어 아크가 사라진다.

## 노드 주파수에는 벌크 간극이 없다 — 아크는 공명이다

phonopy 로 확인: ω = 10.0869 THz 에서 **사영 벌크 연속체가 표면 BZ 창 전체를
덮는다** (커버리지 1.000).  전자계 Weyl 준금속처럼 "간극 속 고립된 아크"가
아니라, 아크가 벌크 연속체와 겹친 **공명**이다.  그래서

- `arc.dat_l` 하나만 rainbow 로 그리면 벌크 배경에 묻혀 아크가 안 보인다.
- 위/아래 표면 차 `dos_l − dos_r` 를 쓰면 벌크가 상쇄되고 아크만 남는다.
  fig 3 (d), fig 4 (c) 가 그것이다.

### 아크임을 확정하는 증거

fig 3 에 겹쳐 그린 두 곡선은 phonopy 로 따로 구한 **band 17 / band 18 의
사영 연속체 경계**다 (각 밴드에 대해 어떤 `k_c` 에서든 ω0 을 지나는 (k_a,k_b)
영역의 경계).  이 두 경계는 **정확히 네 개의 바일 사영점에서 교차**하고,
표면 초과 스펙트럼(빨간 초승달)은 두 경계 사이의 창 안에만 있으며 그 교차점에서
끝난다.  아크가 반드시 만족해야 할 조건이고, 다른 표면 특징은 이걸 만족할 수 없다.

fig 4 는 같은 얘기를 스펙트럼으로 본다.  `k_a = 0.14649` 직선 위에서 벌크 콘
두 개가 `k_b = ±0.07085`, ω = 10.0869 THz 에서 교차한다 — **phonopy 벌크 계산과
Simphony 슬랩 Green 함수가 독립적으로 같은 노드 위치를 준다.**

### 닫힌 고리도 돌려봤다 (fig 로는 안 씀)

`k_a ∈ [0.1015, 0.1915]`, `k_b ∈ [0.0259, 0.1159]` 직사각형(사영점 하나만 포함)
위의 SlabSS 도 계산했다 (`work/slabloop`).  χ=±1 이므로 아크가 반드시 홀수 번
가로질러야 하지만, ω0 에서의 DOS 프로파일은 6.9 → 1.7 → 7.0 으로 매끄럽기만 하고
뾰족한 교차 흔적이 없다.  간극이 없어 아크가 연속체에 묻히기 때문이고,
η = 0.006 THz 로는 분해되지 않는다.  위상적 진술은 맞지만 그림으로 쓰기엔
설득력이 없어서 뺐다.

## ⚠️ Simphony `fermiarc.f90` 의 KPLANE_SLAB 무시 문제

`SlabArc_calc` 은 KPLANE_SLAB 의 두 벡터를 **정수로 올림**한다
(`ceiling(0.6) = 1`).  QPI 가 full BZ 를 필요로 한다는 이유인데, QPI 를 끄고
있어도 항상 그렇게 한다.  그래서 좁은 창을 줘도 실제 계산 영역은
K2D_start 에서 시작하는 **full BZ 한 칸**이다.

이걸 모르고 지정한 범위로 축을 라벨링하면 좌표가 통째로 틀리고, BZ 밖으로
넘어간 부분이 되접혀서 **대칭이 깨진 것처럼 보인다** (실제로 그렇게 한 번
틀렸다).  파일의 kx, ky 폭이 곧 |b1|, |b2| 이므로 그걸로 나누면 reduced 좌표가
나오고, `[-0.5, 0.5)` 로 되감으면 된다.  `fig3_arc.py` 의 `load()` 가 그렇게 한다.
되감은 뒤 대칭 검사: `k_a → −k_a`, `k_b → −k_b` 편차 모두 **0.0000** (범위 5.9).

고치는 패치: `patch/apply_fermiarc_kplane_fix.py` (QPI 일 때만 올림).
이 컨테이너에는 `mpif90` 이 없어 재빌드는 못 했다 — gfortran 문법 검사는
통과했고(직렬 경로 오류 0건), 실제 검증은 클러스터에서 해야 한다.
위 그림들은 **패치 전 바이너리**로 돌린 결과를 full BZ 로 올바르게 해석한 것이라
다시 돌릴 필요는 없다.
