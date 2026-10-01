#!/usr/bin/env python3
"""
cone_fit.py — Weyl 노드 주변 콘을 선형 모형에 맞춰 기울기(type-I/II)와
콘 속도 주축을 뽑는다.

모형:  omega_pm(dk) = omega0 + v . dk  +-  |M dk|        (dk 는 cartesian, 1/Ang)

  - v      : 콘 전체의 기울기 벡터 (두 밴드 평균의 구배)
  - A=M^T M: 콘 계량. 고유값의 제곱근이 세 주축 속도 (THz*Ang)
  - T = sqrt(v^T A^-1 v) : tilt parameter.  T<1 type-I, T>1 type-II

h 를 여러 개 돌려서 T 가 h 에 불변인지 반드시 확인할 것. 노드 좌표가 조금이라도
어긋나 있으면 작은 h 에서 T 가 폭주한다 (그때는 --mode refine 으로 좌표부터 다시 잡는다).

사용법
------
    python3 scripts/cone_fit.py <dir> k1 k2 k3
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph

NOCC = 17          # band 17-18 gap 을 본다 (0-based 16,17)


def main():
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    d = sys.argv[1]
    k0 = np.array([float(x) for x in sys.argv[2:5]])
    ph = get_ph(d, (2, 2, 2))
    B = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi
    Binv = np.linalg.inv(B)
    k0c = k0 @ B

    def fr(kc):
        ph.run_qpoints([list(kc @ Binv)])
        f = ph.qpoints.frequencies[0]
        return f[NOCC - 1], f[NOCC]

    dirs = []
    for a in np.linspace(0, np.pi, 9)[1:-1]:
        for b in np.linspace(0, 2 * np.pi, 13)[:-1]:
            dirs.append([np.sin(a) * np.cos(b), np.sin(a) * np.sin(b), np.cos(a)])
    dirs += [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
    dirs = np.array(dirs)

    for h in [1e-4, 3e-4, 1e-3, 3e-3]:
        rows, rhs, S = [], [], []
        for u in dirs:
            dk = h * u
            w1, w2 = fr(k0c + dk)
            rows.append(dk)
            rhs.append((w1 + w2) / 2)
            S.append((w2 - w1) / 2)
        rows = np.array(rows); rhs = np.array(rhs); S = np.array(S)
        sol, *_ = np.linalg.lstsq(np.hstack([rows, np.ones((len(rows), 1))]), rhs, rcond=None)
        v = sol[:3]
        X = np.array([[p[0]**2, p[1]**2, p[2]**2,
                       2*p[0]*p[1], 2*p[0]*p[2], 2*p[1]*p[2]] for p in rows])
        c, *_ = np.linalg.lstsq(X, S**2, rcond=None)
        A = np.array([[c[0], c[3], c[4]], [c[3], c[1], c[5]], [c[4], c[5], c[2]]])
        ev, evec = np.linalg.eigh(A)
        T = float(np.sqrt(v @ np.linalg.inv(A) @ v))
        print("h=%.0e  T=%.4f (%s)  v_cone=%s" %
              (h, T, "type-II" if T > 1 else "type-I", np.round(np.sqrt(np.abs(ev)), 4)))
        print("        느린 주축 (cartesian) =", np.round(evec[:, 0], 3),
              "   tilt v =", np.round(v, 4))


if __name__ == '__main__':
    main()
