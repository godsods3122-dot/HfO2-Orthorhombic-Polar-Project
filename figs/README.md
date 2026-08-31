# PPT 그림 목록

모든 그림은 300 dpi PNG. 축 규약은 **표준 (a, b, c) — 편극축 = c**로 통일했다
(구조마다 native 축 순서가 달라서 `scripts/figs/axes.py` 가 경로를 매핑한다).
POSCAR 자체는 절대 치환하지 않는다 — FORCE_SETS 와 어긋난다.

| 파일 | 내용 | 만드는 법 |
|---|---|---|
| `fig1_dispersion_parent_pristine.png` | parent_pristine 포논 분산, Γ-X-S-Y-Γ-Z-U-R-T-Z | `run_band.py` → `fig1_2_bands.py` |
| `fig2_bands17_18_node.png` | 밴드 17/18 강조 + 노드 주변 절단 2장 | `fig2_node.py` |
| `fig3_surface_spectrum.png` | ω = 10.0869 THz 표면 스펙트럼 `dos_l − dos_r` (단일 패널, 보조선 없음) | `run_slab.py --mode arc` → `fig3_arc.py` |
| `fig4_slab_arc.png` | `k_b`=0 위의 표면 스펙트럼 (보조) | `run_slab.py --mode ss` → `fig4_slabss.py` |
| `fig5_weyl_cone.png` | 바일 콘 3D + 주축 절단 (type-II 증거) | `fig5_cone.py` |
| `fig6_wcc_4weyl.png` | 4개 노드 각각의 구면 Wilson loop (WCC 합) | Simphony `WeylChirality_calc` → `fig6_wcc.py` |
| `fig7_strain_overlay.png` | biaxial strain 0 / −0.8 / −1.0 / −2.5 / −3.0 % 밴드 중첩 | `run_band.py` ×5 → `fig7_overlay.py` |
| `fig8_weyl_shift.png` | 4점 궤도의 이동 (unstrained → −0.8 %) | `fig8_shift.py` |

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
