# PPT 그림 목록

모든 그림은 300 dpi PNG. 축 규약은 **표준 (a, b, c) — 편극축 = c**로 통일했다
(구조마다 native 축 순서가 달라서 `scripts/figs/axes.py` 가 경로를 매핑한다).
POSCAR 자체는 절대 치환하지 않는다 — FORCE_SETS 와 어긋난다.

| 파일 | 내용 | 만드는 법 |
|---|---|---|
| `fig1_dispersion_parent_pristine.png` | parent_pristine 포논 분산, Γ-X-S-Y-Γ-Z-U-R-T-Z | `run_band.py` → `fig1_2_bands.py` |
| `fig2_bands17_18_node.png` | 밴드 17/18 강조 + 노드 주변 절단 2장 | `fig2_node.py` |
| `fig3_fermi_arc.png` | ω = 10.0869 THz 등주파수 표면 스펙트럼 (bottom / top / bulk) | `run_slab.py --mode arc` → `fig3_arc.py` |
| `fig4_slab_arc_loop.png` | 바일 사영점 하나를 감싸는 닫힌 고리 위의 표면 스펙트럼 | `run_slab.py --mode ss` → `fig4_slabss.py` |
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

## fig 4 경로를 왜 닫힌 고리로 잡았나

직선 절단은 아크를 만난다는 보장이 없다. `k_b = 0` 선도 `k_a = 0` 선도 양쪽
반쪽의 알짜 chirality 가 0 이라 짝수 번(0번 포함) 교차할 수 있다.
반면 **사영점 하나만 감싸는 닫힌 고리**는 χ=±1 이므로 아크가 반드시 홀수 번
가로지른다. 그래서 이 경로가 "아크가 반드시 있는" 경로다.
고리는 `k_a ∈ [0.1015, 0.1915]`, `k_b ∈ [0.0259, 0.1159]` 의 직사각형
(A→B→C→D→A) 이고 안쪽에 `(+0.14649, +0.07085)` 하나만 들어간다.
