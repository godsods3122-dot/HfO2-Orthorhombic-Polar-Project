#!/usr/bin/env python3
"""
sphere_chirality.py — 노드를 둘러싼 구 위에서 Fukui-Hatsugai-Suzuki 격자
Berry flux 를 적분해 chirality 를 구한다. Simphony 와 독립인 검산용.

밴드 1..NOCC 를 한 덩어리로 보고 link variable
    U_ij = det(V_i^dag V_j) / |det|
를 쓴다. 플라케트 flux 는 -arg(U12 U23 U34 U41), 구 전체 합을 2pi 로 나눈 것이
chirality 다.

주의: **전체 부호 규약이 Simphony 와 반대다** (results/node_audit/SUMMARY.md).
절대부호는 Simphony 의 WeylChirality_calc 에 맞추고, 여기서는 상대부호와
"반지름에 대해 정수로 안정한가"만 본다.

반지름을 여러 개 도는 이유: 근처에 위성 노드가 있으면 큰 반지름에서 값이
바뀐다. 실제로 parent_pristine 은 R=0.008 1/Ang 에서 +1 -> -1 로 뒤집힌다
(편극축 방향 dk3=0.0049 위성 두 개가 구 안에 들어온다).

사용법
------
    python3 scripts/sphere_chirality.py <dir> k1 k2 k3
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph

NOCC = 17


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    d = sys.argv[1]
    k0 = np.array([float(x) for x in sys.argv[2:5]])
    ph = get_ph(d, (2, 2, 2))
    B = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi
    Binv = np.linalg.inv(B)
    k0c = k0 @ B

    def vecs(pts):
        ph.run_qpoints([list(p @ Binv) for p in pts], with_eigenvectors=True)
        ev = ph.qpoints.eigenvectors
        return [ev[i][:, :NOCC] for i in range(len(pts))]

    def link(a, b):
        z = np.linalg.det(a.conj().T @ b)
        return z / abs(z)

    nth, nph = 24, 48
    th = np.linspace(0, np.pi, nth + 1)
    phi = np.linspace(0, 2 * np.pi, nph + 1)
    for R in [0.002, 0.004, 0.008]:
        pts = []
        for t in th:
            for p in phi:
                pts.append(k0c + R * np.array([np.sin(t) * np.cos(p),
                                               np.sin(t) * np.sin(p), np.cos(t)]))
        Vl = vecs(pts)
        V = [[Vl[i * (nph + 1) + j] for j in range(nph + 1)] for i in range(nth + 1)]
        tot = 0.0
        for i in range(nth):
            for j in range(nph):
                tot += -np.angle(link(V[i][j], V[i][j + 1]) *
                                 link(V[i][j + 1], V[i + 1][j + 1]) *
                                 link(V[i + 1][j + 1], V[i + 1][j]) *
                                 link(V[i + 1][j], V[i][j]))
        print("R=%.4f 1/Ang   Chern(sphere) = %+.4f" % (R, tot / (2 * np.pi)))


if __name__ == '__main__':
    main()
