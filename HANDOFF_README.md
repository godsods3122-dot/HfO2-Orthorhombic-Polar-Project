# HfO2 위상 포논 — 작업 인수인계

Simphony의 LO-TO 처리 버그를 잡고, 패치된 코드로 HfO2의 Weyl phonon을
조사한 작업의 전체 기록. Claude Code에서 이어서 하기 위한 정리본이다.

---

## 1. 지금까지 확정된 것

### 1.1 Simphony 패치 (`patch/apply_simphony_loto_fix.py`)

원본 Simphony(github.com/fballestermacia/simphony)의 LO-TO 처리에 버그가 있었다.
수정 7건, 검증 완료.

| 검증 항목 | 패치 전 | 패치 후 |
|---|---|---|
| 벌크 밴드 vs phonopy (전 경로 239지점 × 36밴드) | 최대 0.2475 THz | **3.4e-7 THz** |
| Γ에서 acoustic | -0.13 THz | 2e-5 THz |
| slab 허수 모드 | 6개 (최저 -6.43 THz) | **0개** (0.63~24.09 THz) |
| AgP2 공식 예제 회귀 | 1.260234 | 1.260234 (불변) |

수정 내용은 스크립트 상단 docstring에 전부 적혀 있다. 요약하면:

1. **자기항 에르미트화 누락** — QE `rgd_blk`에 있는
   `fmtx = (fmtx + TRANSPOSE(fmtx))/2` 가 이식 과정에서 빠짐.
   결과 dipole 행렬이 비-에르미트(잔차 0.87)가 됨. 벌크/slab 공용 루틴.
2. **q벡터 인덱스 오류** — `q(2)`가 `rec_lattice(2,1)` 대신 `(1,2)`를 참조.
   같은 파일 77~79행에 올바른 판본이 있는 복사-붙여넣기 실수.
   `rec_lattice`가 비대칭인 저대칭 물질에서만 드러나 논문 예제에서는 안 보였음.
3. **phonopy 호환 Gonze 분할** — 새 파일 `gonze_nac.f90`, `LOTO_method` 플래그.
   hr.dat에서 뺀 항과 더하는 항의 Ewald 분할이 같아야 상쇄된다.
   phonopy(c/dynmat.c)와 QE 계보는 위상 규약(G만 vs G+q), Λ 결정 방식,
   dd_q0 대칭화 세 가지가 다르다.
4. **Γ 방향의 k점 순서/MPI 의존성** — `keps`가 암묵적 SAVE라 "직전 k점"을
   재사용. k점이 랭크별 분산이라 `-np N`에 따라 결과가 바뀜.
   phonopy처럼 밴드 경로 구간 방향을 쓰도록 변경.
5. **Γ 비해석항 이중 계산** — gonze가 이미 포함하는 항을 또 더하던 것 제거.
6. **slab 장거리항 이중 계산** — `ek_slab`이 실공간에서 `HmnR += D^LR(R)`
   해두는데 per-k로 또 더함. **원본부터 있던 버그이고 slab 문제의 최대 원인.**
7. **전역 상태 제거** — `added_LR_in_Real_Space`로 공유 `HmnR`을 변형/복원하던
   구조 때문에 한 실행 안에서 계산 순서에 따라 결과가 달라졌음.
   별도 배열 `HmnR_LRslab`으로 분리.

적용:
```bash
cd <simphony 루트>            # src/ 가 보이는 곳
python3 apply_simphony_loto_fix.py
cd src && make clean && make  # module.f90이 바뀌므로 clean 필수
```
`make clean` 생략하면 `readinput.o`가 갱신되지 않아
`Invalid line in namelist PARAMETERS` 오류가 난다.

> 참고: 이 대화에서 쓴 저장소는 `v1.1.0-36-g0d72d28`.
> 사용자 클러스터는 `v1.1.0-34-gf2ae713`이라 `git apply`가 실패하고
> `patch -p1`로 적용됐다. `.rej` 파일이 없었으므로 정상 적용이다.

### 1.2 계산 설정에서 반드시 맞춰야 하는 것

**`SELECTED_OCCUPIED_BANDS` 카드가 `NumOccupied`를 덮어쓴다.**
`input.py` 템플릿은 `1-3`을 쓴다. 17/18 밴드를 보려면 이 카드도
`1-17`로 고쳐야 한다. 이걸 놓쳐서 이 대화에서 chirality가 계속 0으로
나왔다. `NumOccupied`만 바꾸면 소용없다.

`pn.in` 필수 항목:
```
&SYSTEM
  NumOccupied = 17
/
&PARAMETERS
  LOTO_method = 'phonopy'     ! phonopy 유래 hr.dat이면 필수. 기본값은 'qe'
/
SELECTED_OCCUPIED_BANDS
1-17
```

### 1.3 미해결 — Born charge 전치

**어느 쪽이 맞는지 결론 못 냈다.** 상반된 실측이 둘 다 있다.

- 비전치가 맞다는 근거: m2.5에서 Γ→X 경로 X점 주파수가 phonopy 자체 계산과
  비전치일 때 2.808592로 일치, 전치는 2.7939로 어긋남.
- 전치가 맞다는 근거: m1(-1%)에서 전치 `LOTO_BC`로 돌리면 사용자가 보고한
  chirality −1,−1,+1,+1이 정확히 재현됨.

원리적으로는 비전치가 맞아야 한다. phonopy C 코드가 `np.dot(q, born[i])`로
첫 인덱스를 축약하고 논문 Eq.(8)도 같다. 그런데 실측이 반대로 나왔다.

**다만 Weyl의 존재 여부는 전치와 무관하다.** 두 규약 모두 -1%에서 Weyl이
나오고 위치만 다르다:

| 규약 | -1% Weyl 위치 (reduced) | chirality |
|---|---|---|
| 전치 | (±0.08115, 0, ±0.14476) | ∓1, ±1 |
| 비전치 | (±0.06377, 0, ±0.12942) | ∓1, ±1 |

**가려낼 실험**: `SELECTED_OCCUPIED_BANDS 1-17`로 고정한 채 전치만 바꿔
두 번 돌리고, X점 밴드 주파수와 chirality를 **함께** 본다.
둘 다 phonopy와 맞는 조합이 하나만 남아야 정상.

---

## 2. 물리 결과

### 2.1 대칭 구조 (구조마다 축 규약이 다르니 매번 spglib으로 확인할 것)

두 구조 모두 Pca2₁(No.29)이지만 축 순서가 달라 지표 대응이 다르다.
**가정하지 말고 `weyl_scan.py --mode sym`으로 확인할 것.**

| | 거울면 지표 (chirality 0 강제) | 편극축 지표 (Weyl 가능) |
|---|---|---|
| m1 (-1%) | k1, k3 | k2 |
| m2.5 (-2.5%) | k1, k2 | k3 |

이 대화에서 이 대응을 여러 번 착각해서 엉뚱한 평면을 뒤졌다.

### 2.2 -1% (m1): Weyl phonon 존재

비전치+phonopy 조건에서 `(±0.06377, 0, ±0.12942)`, E=10.36596 THz,
gap 5e-15, chirality −1/+1/+1/−1.

`k2=0`(편극축 방향, 거울면 아님)에 있고 `k1`, `k3` 모두 거울면 밖.

### 2.3 -2.5% (m2.5): Weyl phonon 없음

기약 쐐기 전역(41³) + `k3=0` 평면 조밀(161², 간격 0.0031) 스캔 결과,
축퇴는 전부 대칭이 강제하는 자리에만 있다.

- 거울면(`k1` 또는 `k2` = 0/0.5) 위 nodal line — chirality 0
- `k3=0.5` nodal plane **전체** (gap ≤ 2.4e-7, 인접면은 0.061).
  편극축 2₁ screw + 시간역전이 강제.
- **일반 위치 최소 gap = 4.63e-3 THz** (거울면 위는 1e-10 → 7자리 차이)

해상도 확인: -1% Weyl의 `k1`=0.064를 20칸으로 분해하는 격자에서도 골이 없다.
`k2` 방향 미세 스캔(0.001 간격)으로 쌍이 갈라져 있을 가능성도 배제했다.

### 2.4 소멸 메커니즘

소멸하는 쌍은 **c축 수직 glide**로 묶여 있다(m1 기준 `op2 = (k1,k2,−k3)`).
거울이라 chirality를 뒤집으므로 `(k1,0,+k3)`와 `(k1,0,−k3)`는 반대 부호이고,
compressive biaxial strain이 `k3`를 0으로 밀면 glide 불변면에서 만나
**반드시 소멸한다.** 피할 옆길이 없다.

### 2.5 소멸을 피하는 방향 (제안, 미검증)

spglib으로 확인한 결과 **Pca2₁을 보존하는 strain은 수직 성분 3개뿐**이고
어떤 shear도 대칭을 깬다:

| strain | 결과 공간군 |
|---|---|
| 등방 in-plane (εxx=εzz) | Pca2₁ |
| **비등방 in-plane (εxx≠εzz)** | **Pca2₁** |
| c축만 인장 (εzz>0) | Pca2₁ |
| shear εxz | P2₁ |
| shear εxy, εyz | Pc |

따라서 **비등방 in-plane strain**이 유일하게 대칭을 지키면서 쓸 수 있는
새 자유도다. `k3`를 0으로 미는 건 c축 압축이므로, a는 계속 압축하되
c는 덜 누르거나 인장하면 소멸을 늦출 수 있다. 사방정 기판이나 두 방향
격자 부정합이 다른 기판으로 구현 가능.

folding은 권하지 않는다. 단순 재라벨링은 물리를 안 바꾸고, 실제 구조 변조로
주기가 배가되는 경우는 접힌 두 Weyl이 혼성되어 오히려 소멸을 촉진한다.

---

## 3. 다음 할 일

우선순위 순.

1. **Born charge 전치 문제 결론내기** (1.3의 가려낼 실험). 이게 안 정해지면
   모든 좌표값의 신뢰도가 흔들린다.
2. **-2% 구조 계산** — -1%에 있고 -2.5%에 없으니 소멸 지점이 그 사이.
   잔여 gap이 -2.5%에서 매우 작았으므로 소멸 지점이 -2.5% 가까울 수 있다.
   `weyl_scan.py --mode plane --fixed-index <편극축> -N 161` 한 번이면 유무 판정.
3. **비등방 in-plane strain 시험** — εxx는 -1% 고정, εzz를 0%와 +1%로 두 점.
   `k3`가 커지는 방향인지만 보면 전략 타당성이 갈린다.
4. slab 쪽은 패치로 허수 모드가 사라졌지만 외부 참값이 없어 미검증.
   NSLAB을 키웠을 때 내부 층 투영 밴드가 벌크로 수렴하는지가 자체 일관성 검사.

---

## 4. 파일 안내

```
patch/apply_simphony_loto_fix.py   Simphony 패치 (단독 실행, 외부 통신 없음)
scripts/weyl_scan.py               노드 탐색·대칭 분류·정밀화 통합
scripts/extract_weyl_fast.py       Nodes.dat에서 고립점 분리 (8만행 0.6s/85MB)
scripts/phonopy_ref.py             phonopy 참값 밴드 생성
scripts/gonze_ref.py               phonopy Gonze를 numpy로 재현 (이식 검증용)
results/                           밴드 그림, 참값 데이터, 대표 pn.in
```

### 검증용 재현 절차
```bash
# 1) hr.dat + pn.in
python phonopy2TBDAT_LOTO.py POSCAR FORCE_SETS "2 2 2" --born BORN -o wannier90_hr.dat
python input.py --poscar POSCAR --born BORN --hr wannier90_hr.dat   # 전치 문제는 1.3 참조

# 2) pn.in 손보기: LOTO_method='phonopy', NumOccupied=17, SELECTED_OCCUPIED_BANDS 1-17

# 3) 노드 탐색 (phonopy 쪽이 훨씬 빠름)
python scripts/weyl_scan.py --dir . --band 17 --mode sym      # 축 규약 먼저 확인
python scripts/weyl_scan.py --dir . --band 17 --mode wedge -N 41

# 4) 나온 좌표를 pn.in WEYL_CHIRALITY 카드에 넣고 Simphony로 chirality 확인
```

### 주의점
- 노드 탐색은 phonopy로, chirality는 Simphony로. 둘이 3.4e-7 THz 일치하므로
  탐색을 Simphony로 할 이유가 없다.
- `Package = 'Phonopy'` 옵션은 쓰지 말 것. 개발자가 `! UNDER TESTING`이라 명시한
  미검증 분기이고 어떤 릴리스에도 포함된 적이 없다. `LOTO_method`가 우선하므로
  무해하지만 혼란만 준다.
- `LOTO_2D`는 2π 누락을 고쳐 동작은 하나 기본값 `.false.` 유지 권장.
  Sohier 공식은 진공 분리 단일층 전제라 bulk hr.dat을 잘라 만든 slab에 안 맞고,
  실제 slab 문제 원인은 2D가 아니라 이중 계산이었다.
