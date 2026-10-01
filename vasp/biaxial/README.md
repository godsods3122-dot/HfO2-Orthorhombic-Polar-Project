# biaxial strain VASP 입력 일습

```
vasp/biaxial/INCAR_strain      1단계. 변형 축 고정 셀 완화
vasp/biaxial/INCAR_force       2단계. phonopy 유한변위 초격자 힘
vasp/biaxial/INCAR_born        3단계. Born 유효전하 + eps_inf (DFPT)
scripts/make_biaxial_vasp.py   세 단계 디렉터리를 strain 별로 깔아 주는 드라이버
```

## 쓰는 법

```bash
python3 scripts/make_biaxial_vasp.py \
    --ref source/pristine_mirror \
    --strains -3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 \
    --out vasp/biaxial/runs \
    --potcar /path/to/POTCAR
```

각 `runs/eps_*/run.sh` 가 `01_strain → 02_born → 03_force` 를 순서대로 돌고,
끝나면 `out/` 에 **POSCAR, BORN, FORCE_SETS** 가 모인다. 이 셋이
`scripts/weyl_scan.py`, `scripts/phonopy2TBDAT.py` 가 바로 먹는 입력이다.

## 축 규약 — 가정하지 않는다

구조마다 어느 지표가 편극축인지가 다르다 (CLAUDE.md §2.1). 드라이버가
spglib 으로 직접 찾아 **거울축 둘에만 변형을 걸고 편극축은 자유 relax** 로 두며,
`LATTICE_CONSTRAINTS` 줄도 그에 맞춰 써 준다.

| 구조 | 편극축 | 자동 판정되는 LATTICE_CONSTRAINTS |
|---|---|---|
| `pristine_mirror`, `m1_mirror`, `p1_mirror` | a2 | `.FALSE. .TRUE. .FALSE.` |
| `m2.5_mirror`, `m3_mirror`, `parent_pristine`, `hydro_m1` | a3 | `.FALSE. .FALSE. .TRUE.` |

`--lattice-constraints ".TRUE. .FALSE. .TRUE."` 로 덮어쓸 수 있고, 자동 판정과
다르면 어느 축이 고정/자유가 되는지 경고로 풀어서 찍어 준다.

> **VASP 의 의미**: `.TRUE.` = 그 격자벡터를 **완화한다**(자유),
> `.FALSE.` = **고정한다**. 기본값은 `.TRUE. .TRUE. .TRUE.` (= 보통의 ISIF=3).
> 이름이 "CONSTRAINTS" 라 반대로 읽기 쉬우니 주의.
> VASP 6.x 전용 태그다. VASP 5 면 아래 "대체 경로" 참조.

## 반드시 맞춰야 하는 것

- **`GGA`** — 기존 `source/*/POSCAR` 를 만든 범함수와 같아야 한다.
  템플릿은 `PS`(PBEsol)로 두었다. 다르면 strain 계열 전체가 서로 비교 불가다.
- **`ENCUT`** — 세 INCAR 전부 같은 값(600). strain 사이에서도 같아야 한다.
- **`LREAL = .FALSE.`** — 초격자가 커도 `.AUTO.` 로 바꾸지 말 것. 힘이 오염된다.

## 함정

1. **Pulay 응력** — 셀을 움직이는 완화는 한 번에 안 끝난다. `run.sh` 가
   `CONTCAR → POSCAR` 로 2회 재시작하게 되어 있다.
2. **`LEPSILON` 은 `NPAR` 과 함께 못 쓴다.** `INCAR_born` 에 병렬 태그를
   넣지 말 것.
3. **`ISYM = 0` 은 힘 계산에만.** 변위된 초격자를 모체 대칭으로 대칭화하면
   힘상수가 틀어진다. 완화와 Born 에서는 `ISYM = 2` 로 대칭을 지킨다.
4. **`eps_inf` 는 k점 수렴이 느리다.** Born 쪽 메시를 완화보다 조밀하게
   잡는다 (기본 8×8×8 vs 6×6×6).
5. **Born 검산** — `OUTCAR` 의 Z* 를 전 원자에 대해 더하면 0 이어야 한다
   (음향 총합 규칙). `scripts/input.py` 가 이 합을 출력한다.

## 대체 경로 (VASP 5, `LATTICE_CONSTRAINTS` 없음)

`ISIF = 2` 로 이온만 완화한 뒤, 편극축 길이를 몇 점 바꿔 가며
E(c) 를 떠서 최소점을 잡는다. `LATTICE_CONSTRAINTS` 와 같은 결과를 주지만
구조당 5~7회 더 돌아야 한다.
