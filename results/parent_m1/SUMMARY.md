# parent_m1 (신규 데이터셋) — Simphony Weyl phonon 조사

입력: 사용자 제공 `POSCAR(unitcell) + FORCE_SETS + BORN`, `source/parent_m1/`.
patched Simphony 를 이 세션에서 gfortran 으로 빌드해 사용 (`Makefile.gfortran`,
순차 빌드 필수 — `-j` 는 모듈 의존성 때문에 실패한다).

이 문서는 **이 데이터셋만으로** 처음부터 조사한 결과다. 기존 m1/mirror 계열과의
비교는 하지 않았다 (그쪽은 비등방 오염 문제가 있음).

---

## 0. 결론

**band 17-18 에 Weyl phonon 이 없다.**

기약 쐐기 전역 탐색(41³)에서 대칭이 강제하지 않는 축퇴는 **딱 하나**뿐이고,
그것은 선형 Weyl 이 아니라 **2차 접촉(quadratic contact), χ = 0** 이다.

| 위치 (환산) | gap (THz) | 정체 | χ |
|---|---|---|---|
| k3 = 0.5 평면 전체 | ~1e-15 | 편극축 2₁ + 시간역전이 강제하는 **nodal plane** | 0 (강제) |
| k1 또는 k2 = 0, 0.5 | ~1e-10 | 거울면 위 **nodal line** | 0 (강제) |
| **(0.1478247, 0.0723961, 0)** | **5.3e-15** | **2차 접촉** | **0** |

그 접촉점의 모드 주파수는 **10.08787 THz = 41.72 meV = 336.50 cm⁻¹**.

---

## 1. 먼저 — 이 구조는 −1% strain 이 아니다

격자상수 비교 (축 길이 정렬, Å):

| 구조 | | | | 부피 |
|---|---|---|---|---|
| parent_pristine | 4.996952 | 5.018972 | 5.203031 | 130.4897 |
| **parent_m1 (신규)** | **4.997016** | **5.019052** | **5.203462** | **130.5043** |
| m1_mirror (참고: 실제 −1%) | 4.954716 | 5.044581 | 5.167269 | 129.1531 |

**parent_pristine 대비 변형률 = +0.0013%, +0.0016%, +0.0083%.**
즉 사실상 **무변형**이다. 라벨만 −0.01% 로 잘못 적힌 게 아니라 **구조 자체가
−0.01% 수준**이다 (스크립트가 −1% 대신 −0.01% 를 적용한 것으로 보인다).
−1% 결과가 필요하면 다시 만들어야 한다.

## 2. Born charge 전치 문제 — 비전치가 맞다 (측정으로 확정)

`BulkBand_calc` 로 Weyl 근방을 지나는 선(41점)을 뽑아 phonopy 와 전 36밴드 대조:

| LOTO_BC | 전 36밴드 최대 \|Δ\| | band17-18 최소 gap |
|---|---|---|
| 전치 | **1.780e-1 THz** | 1.44e-2 (닫히지 않음) |
| **비전치** | **2.556e-7 THz** | **8.70e-8 (phonopy 9.60e-8 과 일치)** |

레포가 명시한 허용오차 3.4e-7 THz 안에 드는 것은 **비전치뿐**이다.
`scripts/input.py` 의 `--transpose-born` 기본값을 **꺼짐으로 변경**했다.

## 3. 반드시 필요한 pn.in 항목 두 가지

1. **`Package = 'Phonopy'`** — 없으면 `readinput.f90` 기본값이 `'QE'` 이고,
   `readHmnR.f90:160` 이 스피너용 `reorder_wannierbasis` 를 태운다. 고윳값은
   permutation similarity 라 안 변해서 **밴드는 멀쩡해 보이지만** Wannier center
   대응이 깨져 Wilson loop 이 망가진다. `input.py` 가 이 줄을 쓰도록 수정했다.
2. `LOTO_method = 'phonopy'`, `NumOccupied = 17`, `SELECTED_OCCUPIED_BANDS 1-17`
   (기존 CLAUDE.md 지침 그대로).

또 하나: `input.py` 가 만드는 pn.in 에는 **placeholder `WEYL_CHIRALITY` 카드가
이미 들어 있다** (Num_Weyls=2, k=(−0.5,0,0),(0,0.5,0)). 파서는 **첫 번째** 카드를
읽으므로 뒤에 덧붙이면 무시된다. 반드시 **교체**할 것.

## 4. BulkGap_cube_calc 는 이 빌드에서 쓸 수 없다

`BulkGap_cube_calc = T`, `NumOccupied = 17`, KCUBE = [0,0.5]³, 41³ → 17초로 빠르다.
그러나 출력 `GapCube.dat` 의 Ev/Ec 가 **3.54–4.64 THz** 범위다. phonopy(그리고
검증된 `BulkBand_calc`)의 band 17/18 은 **9.8–11 THz** 다.

- `SELECTED_OCCUPIED_BANDS` 카드를 지워도 결과가 완전히 동일 → 그 카드 탓이 아니다.
- Ev/band17 비가 0.357~0.393 로 **일정하지도 않다** → 단순 단위 환산도 아니다.
- 카테시안 열(1–3)을 역격자로 되돌려 같은 k 에서 phonopy 와 대조해도 어떤 고정
  밴드쌍과도 대응되지 않는다. (열 7–9 의 "k1(2pi/a)" 는 환산좌표가 **아니다** —
  이 POSCAR 처럼 축이 치환된 셀에서는 특히 주의.)

**따라서 노드 탐색은 phonopy 로, chirality 만 Simphony 로** 했다 (CLAUDE.md 의
원래 분업 그대로). BulkBand 와 WeylChirality 는 정상 동작을 확인했다.

## 5. 전역 탐색 (phonopy, 기약 쐐기 41³)

```
쐐기 최소 gap = 1.84e-10
gap<0.02 국소최소 35개 → 거울면 위/수렴 30개, 일반위치 5개
```

일반위치 후보 5개를 Simphony `WeylChirality_calc` (비전치, 반지름 0.002) 로 확인:

| k (환산) | Simphony χ | 판정 |
|---|---|---|
| (0.4230957, 0.0919059, **0.5**) | −1 | k3=0.5 **nodal plane 위** → 구가 축퇴면을 관통, **인공물** |
| (0.3201510, 0.3048336, **0.5**) | +1 | 동일하게 인공물 |
| (0.2319125, 0.3747003, **0.5**) | +1 | 동일하게 인공물 |
| **(0.1478247, 0.0723961, 0.0)** | **0** | 유일한 진짜 일반위치 |
| (0.0121886, 0.4478596, **0.5**) | 0 | nodal plane 위 |

k3=0.5 는 편극축 2₁ screw + 시간역전이 band17-18 축퇴를 **평면 전체에** 강제하는
자리다 (쐐기 스캔에서 그 평면 위 여러 점의 gap 이 모두 1e-15~1e-10). 축퇴면을
관통하는 구에서 계산한 Wilson loop 은 정의되지 않으므로 거기 나온 ±1 은 의미가 없다.

## 6. 유일한 일반위치 축퇴는 Weyl 이 아니다

`(0.1478247, 0.0723961, 0)` 주변 구면에서 밴드간 최소 gap:

| R (환산) | 16-17 | **17-18** | 18-19 |
|---|---|---|---|
| 0.001 | 8.82e-1 | **1.347e-5** | 3.07e-1 |
| 0.002 | 8.77e-1 | **5.354e-5** | 3.05e-1 |
| 0.003 | 8.72e-1 | **1.203e-4** | 3.04e-1 |
| 0.004 | 8.67e-1 | **2.139e-4** | 3.02e-1 |
| 0.008 | 8.47e-1 | **8.560e-4** | 2.95e-1 |

- 아래위 밴드(16-17, 18-19)는 0.88 / 0.30 THz 로 **충분히 떨어져 있다** →
  17번 밴드만의 Berry flux 도 원리적으로는 정의된다.
- 그런데 17-18 gap 이 **정확히 R² 로 자란다** (1 : 3.97 : 8.9 : 15.9 ≈ 1:4:9:16).
  선형 Weyl 이면 R 에 비례해야 한다. → **2차 접촉**이다.
- 절대 크기도 극히 작다: R=0.008 에서도 8.6e-4 THz.

이 때문에 단일 밴드 Berry flux 는 반지름에 따라 부호가 뒤집힌다
(R≤0.002 에서 −1, R≥0.004 에서 +1, 격자 24×48 과 48×96 에서 동일 → 격자
수렴 문제가 아니라 gap 이 1e-5 THz 라 고유벡터가 수치적으로 불안정한 것).
±0.007 상자를 25³ 로 훑어도 그 안에 노드는 **이것 하나뿐**이라, 근처에 숨은
반대부호 짝이 있어서 생긴 현상도 아니다.

**Simphony 의 17밴드 Wilson loop 은 χ = 0 을 준다.** 2차 분산 + 다중밴드
Wilson loop 0 → **단극 전하 0**. Weyl 이 아니다.

> 신뢰도 단서: gap 이 1e-5 THz 수준이라 이 판정은 수치적으로 민감하다.
> 확정하려면 더 조밀한 FORCE_SETS(3×3×3 supercell)나 더 타이트한 DFT 수렴이 낫다.

## 7. 다음

1. **진짜 −1% 구조를 다시 만들어 재계산** (2절). 지금 것은 무변형이다.
2. 무변형 parent 에 Weyl 이 없다는 것 자체는 정보다 — 변형을 걸어야 나온다는 뜻.
   부호와 크기를 정하려면 −1%, −2% 두 점이면 충분하다.
3. 2차 접촉이 변형에 따라 두 개의 선형 Weyl 로 갈라지는지 보는 것이 가장 자연스러운
   다음 질문이다 (2차 접촉은 χ=±2 double Weyl 또는 χ=0 이 갈라지기 직전 상태다).

## 재현

```bash
cd simphony/src && cp Makefile.gfortran Makefile && make clean && make   # -j 금지
python3 scripts/phonopy2TBDAT.py POSCAR FORCE_SETS "2 2 2" --born BORN -o wannier90_hr.dat
python3 scripts/input.py --poscar POSCAR --born BORN --hr wannier90_hr.dat   # 비전치가 기본값
# pn.in: LOTO_method='phonopy', NumOccupied=17, SELECTED_OCCUPIED_BANDS 1-17,
#        WeylChirality_calc=T, 그리고 placeholder WEYL_CHIRALITY 카드를 '교체'
python3 scripts/weyl_scan.py --dir source/parent_m1 --band 17 --mode wedge -N 41
simphony/src/pn.x
```
