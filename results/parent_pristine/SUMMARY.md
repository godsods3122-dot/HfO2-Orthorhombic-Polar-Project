# parent_pristine — Weyl phonon 있음 (Simphony 확정)

Simphony(패치본, 이 세션 gfortran 빌드)로만 조사. phonopy 는 노드 탐색과 교차검증에만 사용.
외부 문헌 참조 없음. `source/parent_pristine/`.

---

## 0. 결론

**band 17-18 에 진짜 Weyl phonon 이 있다.**

`(k1, k2) ≈ (0.1464, 0.0709)` 에 **k3 방향으로 3개가 일렬로** 늘어서 있다.
모드 주파수는 셋 다 **10.087 THz = 41.72 meV = 336.5 cm⁻¹**.

| k (환산) | gap (THz) | χ (Simphony) | χ (Berry flux) |
|---|---|---|---|
| (0.1464400, 0.0710381, **+0.0049255**) | 0 | **−1** | +1 |
| (0.1464927, 0.0708493, **0**) | 0 | (0, 아래 3절) | −1 |
| (0.1464400, 0.0710381, **−0.0049255**) | 0 | **−1** | +1 |
| **세 개를 감싸는 구 (R=0.010)** | — | **−1** | **+1** |

두 방법의 전체 부호 규약이 반대다(크기·상대부호는 일치). **합 규칙이 양쪽에서 맞는다**:
−1 + (+1) + (−1) = −1 (Simphony), +1 + (−1) + (+1) = +1 (Berry flux).

거울지표가 k1, k2 이므로 각 노드는 **4점 궤도**를 이루고 χ 가 교대한다.
Simphony 직접 확인 (Nk=100, R=0.002):

```
( 0.14644,  0.07104, 0.00493) -> -1
( 0.14644, -0.07104, 0.00493) -> +1
(-0.14644,  0.07104, 0.00493) -> +1
(-0.14644, -0.07104, 0.00493) -> -1
( 0.14644,  0.07104,-0.00493) -> -1
```

즉 이 계열은 BZ 전체에 **12개**(3 × 4점 궤도)다.

## 1. 해밀토니안 검증 — Simphony = phonopy

세 노드를 관통하는 선(k3 = −0.01 → +0.01, 81점)에서 `BulkBand_calc` 대 phonopy:

```
전 36밴드 최대 |simphony − phonopy| = 2.294e-07 THz
band17-18 국소최소  simphony: k3=±0.005 에서 6.83e-06
                   phonopy : k3=±0.005 에서 6.83e-06   (동일)
```

레포 허용오차 3.4e-7 THz 안에 든다. **비전치 LOTO_BC + `Package = 'Phonopy'`** 조합이며,
`scripts/input.py` 가 이제 이 둘을 기본으로 낸다.

## 2. 전역 탐색 (기약 쐐기 41³)

```
쐐기 최소 gap = 2.01e-10
gap<0.02 국소최소 36개 → 거울면 위/수렴 33개, 일반위치 3개
   (0.1464400, 0.0710381, 0.0049255)  gap=8.9e-15   <- 위 Weyl 계열
   (0.2499400, 0.1370633, 0.5000000)  gap=8.8e-09   <- k3=0.5 nodal plane
   (0.4683408, 0.4643373, 0.5000000)  gap=1.8e-08   <- k3=0.5 nodal plane
```

k3 = 0.5 는 편극축 2₁ screw + 시간역전이 축퇴를 평면 전체에 강제하는 자리다.
Simphony 도 그 둘에 대해 **χ = 0** 을 준다 (강제된 값과 일치).

**따라서 band 17-18 의 Weyl 은 위 한 계열뿐이다.**

## 3. ⚠️ 중요한 함정 — `Nk1`, `Nk2` 를 100 이상으로 둘 것

이 노드들은 매우 **무르다**. 반지름 0.0015 인 구 위에서 gap 이 겨우 ~1e-4 THz 다.
Simphony 의 구면 Wilson loop 은 `Nk1 × Nk2` 로 이산화되는데, 격자가 성기면
감김을 놓치고 **조용히 0 을 돌려준다**:

| Nk1 = Nk2 | 결과 |
|---|---|
| 15 | 전부 0 ← **거짓 음성** |
| 50 | 전부 0 ← **거짓 음성** |
| **100** | **−1, +1, +1, −1, −1** (정상) |

`input.py` 기본값이 50 이므로 **그대로 쓰면 Weyl 을 놓친다.**

바깥 두 노드는 Nk=100 에서 풀리지만, **가운데(k3=0) 노드는 Nk=200, R=0.001 에서도
Simphony 가 0 을 준다.** 그런데 세 개를 감싸는 구(R=0.010, Nk=200)에서는 −1 이 나온다.
가운데가 0 이라면 합이 −2 여야 하므로, **Simphony 자신의 합 규칙이 가운데 노드가
0 이 아님을 말해 준다.** 독립 Berry flux 도 반지름 0.0005~0.002, 격자 32×64 / 48×96
전부에서 −1 로 흔들림 없다. 가운데 노드의 개별 값 0 은 Simphony 쪽 수렴 실패다.

## 4. 이전 세션 결과가 왜 틀렸었나

`results/weyl_trend/SUMMARY.md` 는 parent_pristine 에 대해 "후보 3개, chirality 전부 0"
이라고 적었다. 그 계산은 **전치 LOTO_BC** 였고 `Package` 줄도 없었다. 이번에 측정으로
확인된 바:

- 전치 LOTO_BC 는 전 36밴드에서 phonopy 와 1.8e-1 THz 어긋나고 Weyl 자리에서 gap 이
  닫히지도 않는다 (비전치는 2.3e-7 일치).
- `Package` 줄이 없으면 기본값 `'QE'` 라 `readHmnR` 이 스피너용 `reorder_wannierbasis`
  를 태운다. 고윳값은 permutation 이라 안 변해 **밴드는 멀쩡해 보이지만** Wilson loop
  이 망가진다.
- 거기에 `Nk1=Nk2=50` 의 거짓 음성까지 겹쳤다.

세 가지가 전부 **"조용히 0"** 으로 수렴하는 실패 양상이라 알아채기 어려웠다.

## 5. 발표용 한 줄

> 무변형 parent Pca2₁ HfO₂ 의 band 17-18 에 **χ = ±1 Weyl phonon 12개**가 존재한다.
> (k1,k2) ≈ (0.146, 0.071) 에서 편극축 방향으로 세 개가 Δk3 = 0.0049 간격으로 늘어서
> 있고 알짜 전하는 ±1, 에너지는 **41.72 meV (336.5 cm⁻¹)**.

## 재현

```bash
cd simphony/src && cp Makefile.gfortran Makefile && make clean && make    # -j 금지
cd work && python3 scripts/phonopy2TBDAT.py POSCAR FORCE_SETS "2 2 2" --born BORN -o wannier90_hr.dat
python3 scripts/input.py --poscar POSCAR --born BORN --hr wannier90_hr.dat
# pn.in: WeylChirality_calc=T, NumOccupied=17, SELECTED_OCCUPIED_BANDS 1-17,
#        LOTO_method='phonopy',  Nk1=Nk2=100 이상,
#        placeholder WEYL_CHIRALITY 카드를 '교체'(뒤에 덧붙이면 무시됨)
simphony/src/pn.x
```
