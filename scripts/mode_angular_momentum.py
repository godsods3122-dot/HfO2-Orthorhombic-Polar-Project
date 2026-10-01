#!/usr/bin/env python3
"""
mode_angular_momentum.py — 포논 모드의 원편광(타원)과 각운동량.

"O 가 원을 그리고 그 평면이 도는 것처럼 보인다" 를 정량화한다.

원리
----
질량가중 고유벡터를 e_j = a_j + i b_j 라 하면 원자 j 의 실제 변위는

    u_j(t) = Re[ (e_j / sqrt(m_j)) e^{-i w t} ]

이고, 이는 **a_j 와 b_j 가 張하는 고정된 평면** 안의 타원을 그린다.
한 모드·한 q 에서 타원면은 시간에 따라 돌지 않는다 — 이건 수학적으로 고정이다.
시간평균 각운동량은

    <l_j> = 2 (a_j x b_j)      [hbar 단위, 질량가중 기준]

원형도  circ = 2|a x b| / (|a|^2 + |b|^2) 는 1 이면 완전한 원, 0 이면 직선 진동.

"평면이 도는" 인상의 출처는 셋 중 하나다:
  (1) 원자마다 타원면이 다르다 (HfO2 에서는 법선이 최대 ~89도까지 벌어진다),
  (2) 여러 단위포를 동시에 그리면 q.R 위상 때문에 실공간에서 나선이 된다
      — 원편광 포논이 진행할 때의 진짜 모습이다,
  (3) 축퇴 근처에서 두 모드를 섞으면 맥놀이로 천천히 돈다
      (단, 노드 근처 Δw ~ 0.003 THz vs w ~ 10 THz 라 진동 3000여 주기에 한 번이다).

사용법
------
    python3 scripts/mode_angular_momentum.py <dir> k1 k2 k3 [mode0 mode1 ...]

mode 는 0-based (band 17 을 보려면 16). 생략하면 16, 17.
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    d = sys.argv[1]
    k = [float(x) for x in sys.argv[2:5]]
    modes = [int(x) for x in sys.argv[5:]] or [16, 17]

    ph = get_ph(d, (2, 2, 2))
    sym = ph.primitive.symbols
    ph.run_qpoints([k], with_eigenvectors=True)
    f = ph.qpoints.frequencies[0]
    ev = ph.qpoints.eigenvectors[0]

    print("%s  k=(%.6f, %.6f, %.6f)" % (d, *k))
    for m in modes:
        e = ev[:, m].reshape(-1, 3)
        a, b = e.real, e.imag
        L = 2 * np.cross(a, b)
        print("\n  mode %d (band %d)  w = %.6f THz   총 각운동량 L = (%+.4f, %+.4f, %+.4f) hbar"
              % (m, m + 1, f[m], *L.sum(axis=0)))
        print("    원자   |a|     |b|    원형도   타원면 법선(cart)          |l_j|")
        ns, ws = [], []
        for j in range(len(e)):
            na, nb = np.linalg.norm(a[j]), np.linalg.norm(b[j])
            cr = np.cross(a[j], b[j])
            ncr = np.linalg.norm(cr)
            den = na**2 + nb**2
            circ = 2 * ncr / den if den > 1e-12 else 0.0
            n = cr / ncr if ncr > 1e-12 else np.zeros(3)
            if den > 1e-4:
                ns.append(n); ws.append(den)
            print("    %-4s %.4f  %.4f   %.3f    (%+.3f,%+.3f,%+.3f)    %.4f"
                  % (sym[j], na, nb, circ, *n, np.linalg.norm(L[j])))
        if len(ns) > 1:
            ns = np.array(ns)
            ref = ns[int(np.argmax(ws))]
            ang = np.degrees(np.arccos(np.clip(np.abs(ns @ ref), 0, 1)))
            print("    타원면 법선 흩어짐: 최대 %.1f도, 평균 %.1f도 (진폭 가장 큰 원자 기준)"
                  % (ang.max(), ang.mean()))


if __name__ == '__main__':
    main()
