# PPT 그림 목록

모든 그림은 300 dpi PNG. 축 규약은 **표준 (a, b, c) — 편극축 = c**로 통일했다
(구조마다 native 축 순서가 달라서 `scripts/figs/axes.py` 가 경로를 매핑한다).
POSCAR 자체는 절대 치환하지 않는다 — FORCE_SETS 와 어긋난다.

| 파일 | 내용 | 만드는 법 |
|---|---|---|
| `fig1_dispersion_parent_pristine.png` | parent_pristine 포논 분산, Γ-X-S-Y-Γ-Z-U-R-T-Z | `run_band.py` → `fig1_2_bands.py` |
| `fig2_bands17_18_node.png` | 밴드 17/18 강조 + 노드 주변 절단 2장. 안쪽 Γ 는 LO-TO 방향 의존으로 불연속이라 선을 끊어 그린다 (아래 주의). Γ-X 구간에 **Weyl 노드 에너지**(10.087 THz)를 반투명 띠로, 그 구간에서 밴드 17/18 이 실제로 만나는 점을 원으로 표시 | `fig2_node.py` |
| `fig3_surface_spectrum.png` | ω = 10.0869 THz 표면 스펙트럼 `dos_l − dos_r` (단일 패널, 보조선 없음) | `run_slab.py --mode arc` → `fig3_arc.py` |
| `fig4_slab_arc.png` | `k_b`=0 위의 표면 스펙트럼 (보조) | `run_slab.py --mode ss` → `fig4_slabss.py` |
| `fig5_weyl_cone.png` | 바일 콘 3D + 주축 절단 (type-II 증거) | `fig5_cone.py` |
| `fig6_wcc_4weyl.png` | 4개 노드 각각의 구면 Wilson loop (WCC 합) | Simphony `WeylChirality_calc` → `fig6_wcc.py` |
| `fig7_strain_overlay.png` | biaxial strain 0 / −0.8 / −1.0 / −2.5 / −3.0 % 밴드 중첩 | `run_band.py` ×5 → `fig7_overlay.py` |
| `fig8_weyl_shift.png` | 4점 궤도의 이동 (unstrained → −0.8 %) | `fig8_shift.py` |

## 발표(PPT)용 단일 패널 — 폰트 크게, 패널 하나씩

| 파일 | 내용 | 만드는 법 |
|---|---|---|
| `fig2_path_ppt.png` | fig 2 의 (a) 표준 경로 밴드만. 고대칭점 문자 60 pt (볼드 없음), y 레이블 48 pt | `fig2_path_ppt.py` |
| `fig2_cut_ka.png` | 노드를 지나는 **$k_a$ 방향** 절단. 축은 Δk 가 아니라 절대 $k_a$ (노드 0.14649 를 지난다). 교점 레이블은 `χ = +1` 만 | `fig2_cut_ppt.py` |
| `fig2_cut_kb.png` | 같은 것의 **$k_b$ 방향** 절단 (절대 $k_b$, 노드 0.07085) | 〃 |
| `fig5a_weyl_cone.png` | fig 5 의 (a) 3D 콘만, 폰트 크게, 레이블 `χ = +1` 만 | `fig5_cone.py` (같이 나온다) |

fig2 의 (b)=$k_a$ 절단, (c)=$k_b$ 절단이다. 발표에는 X자로 갈라지는 $k_b$ 쪽이
읽기 쉽고, type-II(두 가지가 같이 올라감)를 보이려면 $k_a$ 쪽이 맞다.

## 숫자 요약 (발표에 그대로 쓸 수 있는 값)

| 항목 | 값 | 출처 |
|---|---|---|
| 바일 4점 위치 (k_c=0 평면) | `(±0.14649, ±0.07085, 0)` | Simphony BulkGap 401² + 정밀화 |
| 주파수 | 10.0869 THz = 41.72 meV | 〃 |
| chirality | `+1, −1, −1, +1` | Simphony `WeylChirality_calc` (PN.out) |
| gap (노드에서) | 6.5e-8 THz | phonopy |
| **Γ-X 위 교점** | `t=0.185742` (Γ→X 의 37.15 %), **10.192507 THz**, gap 1.6e-11 | phonopy 정밀화 |
| tilt 파라미터 T | 2.748 → **type-II** | 3D 국소 전개 |
| 면내 콘 기울기 | 1.1192 / 0.7619 THz/rlu | 〃 |
| 편극축(c) 방향 속도 | ≈0.0032 THz/rlu (면내의 1/350) | 〃 |
| −0.8 % 에서 위치 | `(±0.12942, ±0.06377, 0)` | 〃 |
| 이동량 | `Δk_a = −0.0171, Δk_b = −0.0071, |Δk| = 0.0185` | 〃 |
| 소멸 구간 | −0.8 % 와 −1.0 % 사이 | BulkGap + phonopy 둘 다 |

## ⚠️ fig 2 — 안쪽 Γ 에서 밴드를 이어 그리면 안 된다

LO-TO 비해석항은 Γ 로 **접근하는 방향**에 의존하므로 `Y→Γ` 의 Γ 값과 `Γ→Z` 의 Γ 값이
다르다. 물리이지 수치 오차가 아니다.

| 밴드 | 좌 | 우 | 점프 |
|---|---|---|---|
| band 20 | 10.2335 | 11.5486 | **1.3151 THz** |
| band 36 | 23.5481 | 22.0976 | 1.4506 |
| band 17 | 9.7933 | 9.8198 | 0.0265 |
| band 18 | 10.0198 | 9.9835 | 0.0363 |

0.1 THz 넘게 튀는 밴드가 **13개**다. 이어 그리면 거의 수직인 가짜 선이 생긴다
(특히 band 20 은 그림 세로 범위 한가운데를 가로지른다).

**끊는 방법**: 그 인덱스에 NaN 하나를 끼운다. 값 자체가 불연속이라 그것으로 충분하다.
끊는 곳은 **안쪽 Γ 한 군데뿐**이다. (틈을 인위적으로 벌리는 판본도 만들어 봤으나
축이 왜곡되어 되돌렸다.)

## fig 2 의 Γ-X 구간 표시 두 가지 — 서로 다른 것이다

**반투명 붉은 띠 = Weyl 노드 에너지 10.086936 THz.** 노드 자체는 `(0.1465, 0.0708, 0)`
로 **경로 위에 없으므로** 점으로 찍을 수 없다. 그 대신 "노드가 어느 높이에 있는가" 를
Γ-X 구간에 띠로 표시한다.

**빨간 원 = 밴드 17/18 이 Γ-X 위에서 실제로 만나는 점.** 이건 노드와 다른 것이고
높이도 다르다 (10.193 vs 10.087 THz). Γ-X 는 거울면 위 선이라 χ=0 이 강제되므로
여기 만나는 점은 Weyl 이 될 수 없다.
**t = 0.185742, E = 10.192507 THz, gap 1.6e-11**.

격자 데이터 `bulkek_parent_pristine.dat` (Γ-X 를 120점) 의 최소 gap 은 1.9e-3 THz 인데,
이건 격자가 교점을 비껴간 것이다. 그림의 표시 위치는 phonopy 정밀화 값을 쓴다.

참고로 이 경로 위에는 **선 전체가 축퇴인 구간**도 있다 — `X-S` 전 구간과
`Z-U-R-T-Z` 전 구간 (후자가 편극축 nodal plane). 점이 아니므로 마커로 표시하지 않았다.

## 슬랩 설정 (fig 3, 4)

```
SURFACE  (1 0 0) / (0 1 0)     → 표면 법선 ∥ c = 편극축
NP = 2                          hr.dat 의 R3 범위가 ±1 이라 2 면 정확하다
LOTO_method = 'phonopy'         (빠뜨리면 밴드가 0.03 THz 조용히 어긋난다)
NumOccupied = 17,  SELECTED_OCCUPIED_BANDS 1-17
```

표면 법선을 c 로 잡아야 4개 사영점이 서로 겹치지 않는다. a 나 b 를 법선으로
잡으면 부호가 반대인 두 점이 같은 자리에 사영되어 아크가 사라진다.

## 페르미 아크 — Simphony 버그 두 개를 고치고 나서야 보인다

`patch/` 에 있는 세 스크립트를 다 적용해야 fig 3 이 제대로 나온다.
자세한 내용은 `results/parent_pristine/SUMMARY.md` 17절.

| 패치 | 무엇 |
|---|---|
| `apply_fermiarc_loto_fix.py` | **SlabArc 가 LO-TO 를 건너뛴다.** `fermiarc.f90:214` 가 조건 없이 plain `ham_qlayer2qlayer` 를 부른다 (`surfstat.f90:142` 는 제대로 갈래를 탄다). fig 3 만 해당, fig 4 는 원래 정상 |
| `apply_fermiarc_kplane_fix.py` | KPLANE_SLAB 벡터를 정수로 올려(ceiling) 항상 full BZ 를 계산한다 |
| `apply_surfdos_only_norm_fix.py` | `dos_l_only` 가 항상 eps9. `dos_l` 은 궤도 36개, `dos_bulk` 은 72개를 합하는데 그냥 뺀다 |

### 왜 raw DOS 로는 안 보이나

노드 주파수에서 **사영 벌크 연속체가 표면 BZ 전체를 덮는다** (phonopy 로 확인,
커버리지 1.000).  전자계 Weyl 준금속처럼 "간극 속 고립된 아크"가 아니라
벌크와 겹친 **공명**이다.  두 가지 방법으로 벌크를 걷어낸다:

- **표면 전용** `ρ_surf − ρ_bulk·N_top/N_dim` — fig 4 (c)
- **위/아래 표면 차** `dos_l − dos_r` — 벌크가 정확히 상쇄된다. fig 3 (d)

### ⚠️ 식별 가능한 아크는 없다 (이전 주장 철회)

`dos_l − dos_r` 의 밝은 능선(`k_a ≈ ±0.10`)을 한때 아크로 봤으나 아니다.
능선이 노드의 `k_b = ±0.07085` 를 그냥 통과해 `|k_b| ≈ 0.11` 까지 이어지고,
그 위치는 band 17 과 18 이 **둘 다** 사영된 영역이다.  아크는 사영점에서
끝나는 열린 곡선이어야 하고 한 밴드만 사영된 영역에 있어야 한다.
실제 아크 영역의 표면 무게는 배경과 구분되지 않는다 (중앙값 0.444 vs 0.459).

원인은 위와 같다: 노드 주파수에 벌크 간극이 없어 아크가 있더라도 공명이고
현재 분해능으로 분리되지 않는다.  자세한 근거는 SUMMARY 17.5.

## ⚠️ 미해결 — 슬랩 LO-TO 가 표면 거울 대칭을 깬다

이 표면은 z 이동이 없는 a-glide(⊥b)를 그대로 가지므로
`dos(k_a, k_b) = dos(k_a, −k_b)` 여야 한다.  그런데 반전(시간역전)은
`0.000e+00` 로 정확한 반면 거울은 log DOS 중앙값 0.024 / 최대 0.66
(전체 범위 7.8) 어긋난다.  벌크 phonopy 는 같은 거울을 1e-7 THz 로 지키고,
SlabSS 도 같은 크기로 어긋나므로 위 패치들이 만든 것이 아니다.
`ham_qlayer2qlayer_LOTO` 쪽으로 보이나 **원인은 아직 못 짚었다.**
그림에서는 알려진 대칭으로 평균만 냈다.

## 빌드 참고

이 컨테이너에는 `mpif90` 이 없어 `src/Makefile.gfortran` 으로 직렬 빌드했다.
검증: 같은 SlabSS 입력에서 기존 MPI 바이너리와 `dos.dat_l/_r/_bulk` 모두
`max|Δ| = 0.000e+00` 으로 완전 일치.
