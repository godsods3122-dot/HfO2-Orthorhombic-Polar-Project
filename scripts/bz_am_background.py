#!/usr/bin/env python3
"""
bz_am_background.py — BZ 전체의 모드 각운동량/자기모멘트 배경.

"노드 하나"가 아니라 **bulk 배경**을 본다. 묻는 것은 셋이다:
  1. 관심 에너지 근처에 몇 개의 가지가 깔려 있는가 (에너지 분해 측정의 배경),
  2. 원편광(|L|)이 노드에서만 특별한가, 아니면 결정 전체에 흔한가,
  3. 그 배경이 알짜로 남는가 (시간반전이 지우는가).

각운동량  L = 2 sum_j (a_j x b_j)  [hbar, 질량가중 고유벡터]
자기모멘트 mu_a/mu_B = -sum_j (m_e/m_j) eps_abc Z*_j,cd Im[e_jb e_jd^*]
           (강체이온 + Born 유효전하. 전자 매개 기여는 포함하지 않는다.)

격자는 Γ 를 피한 Monkhorst-Pack 이고 ±q 짝이 모두 들어 있으므로,
시간반전 L(-q) = -L(q) 가 성립하면 알짜 합은 정확히 0 이어야 한다.
합/개별합 비가 그 검산이다.

사용법
------
    python3 scripts/bz_am_background.py <dir> <N> <E0_THz>

    python3 scripts/bz_am_background.py source/parent_pristine 20 10.0868
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph

ME_U = 5.48579909e-4       # 전자질량 [u]


def main():
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    d = sys.argv[1]
    N = int(sys.argv[2])
    E0 = float(sys.argv[3])

    ph = get_ph(d, (2, 2, 2))
    m = np.array(ph.primitive.masses)
    Z = np.array(ph.nac_params['born'])
    nat = len(m)
    nb = 3 * nat
    eps = np.zeros((3, 3, 3))
    for a, b, c in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        eps[a, b, c] = 1; eps[a, c, b] = -1

    g = (np.arange(N) + 0.5) / N - 0.5
    QS = np.array([[x, y, z] for x in g for y in g for z in g])
    print("%s   BZ 표본 %d점 x %d모드" % (d, len(QS), nb))

    F = np.zeros((len(QS), nb))
    LV = np.zeros((len(QS), nb, 3))
    MU = np.zeros((len(QS), nb, 3))
    CH = 400
    for s in range(0, len(QS), CH):
        q = QS[s:s+CH]
        ph.run_qpoints([list(x) for x in q], with_eigenvectors=True)
        F[s:s+CH] = ph.qpoints.frequencies
        e = ph.qpoints.eigenvectors.transpose(0, 2, 1).reshape(len(q), nb, nat, 3)
        LV[s:s+CH] = 2 * np.cross(e.real, e.imag).sum(axis=2)
        G = np.einsum('qmjb,qmjd->qmjbd', e, e.conj()).imag
        MU[s:s+CH] = -np.einsum('j,abc,jcd,qmjbd->qma', ME_U/m, eps, Z, G)
        if s % 2000 == 0:
            print("  %d/%d" % (s, len(QS)), flush=True)

    Lm = np.linalg.norm(LV, axis=2)
    Mm = np.linalg.norm(MU, axis=2)

    print("\n=== 1. %.4f THz 주변에 깔린 가지 ===" % E0)
    for w in (0.05, 0.10, 0.25, 0.50):
        sel = np.abs(F - E0) < w
        bands = [int(b) + 1 for b in sorted(set(np.where(sel)[1]))]
        print("  ±%.2f THz: 모드 %6d개 (전체의 %.2f%%), 관여 밴드 %s, q점의 %.1f%%"
              % (w, sel.sum(), 100*sel.mean(), bands, 100*sel.any(axis=1).mean()))

    print("\n=== 2. |L| 분포 (hbar) ===")
    def stat(x, lbl):
        print("  %-24s 중앙값 %.3f  평균 %.3f  90%%분위 %.3f  최대 %.3f"
              % (lbl, np.median(x), x.mean(), np.percentile(x, 90), x.max()))
    stat(Lm.ravel(), "전 밴드 x 전 BZ")
    sel = np.abs(F - E0) < 0.10
    if sel.any():
        stat(Lm[sel], "E0 ±0.10 THz")

    print("\n=== 3. 시간반전 상쇄 ===")
    print("  Σ L  = (%+.3e, %+.3e, %+.3e) hbar   Σ|L|  = %.1f" % (*LV.sum(axis=(0, 1)), Lm.sum()))
    print("  Σ mu = (%+.3e, %+.3e, %+.3e) mu_B   Σ|mu| = %.3e" % (*MU.sum(axis=(0, 1)), Mm.sum()))
    print("\n=== 4. 가지 수 ===")
    print("  원자 %d x 3 = %d,  phonopy 가 준 가지 수 = %d  -> 여분의 가지 없음"
          % (nat, 3*nat, F.shape[1]))


if __name__ == '__main__':
    main()
