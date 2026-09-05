#!/usr/bin/env python3
"""
arc_projection.py — Weyl 궤도의 chirality 를 구하고, 어느 표면/계면 배향에서
포논 arc 가 나오는지 판정한다.

원리
----
법선 n 인 표면/계면의 2D BZ 로 Weyl 들을 투영했을 때, 투영점의 **알짜 chirality**
가 0 이 아닌 곳에서만 arc 가 시작·끝난다. 반대 부호 Weyl 이 같은 점에 겹쳐
투영되면 알짜가 0 이라 arc 가 없다.

결과 (pristine_mirror, band17-18, 편극축 = b)
-------------------------------------------
  법선 ∥ a  → 알짜 0  → arc 없음
  법선 ∥ c  → 알짜 0  → arc 없음
  법선 ∥ b (편극축) → 네 투영점 모두 알짜 ±1 → **arc 있음**

즉 이 Weyl 계열의 arc 는 **편극축에 수직인 계면에서만** 보인다. 그건 강유전
커패시터의 **전극 계면**이다 — 소자가 이미 갖고 있는 계면이고, 대칭을 아무것도
바꾸지 않아도 된다.

180° 도메인 벽은 전기적으로 중성이려면 법선이 a 또는 c 여야 하는데, 그 둘은
알짜 0 이라 **도메인 벽 자체에는 이 계열의 arc 가 없다.**

또한 두 180° 도메인은 모든 k 에서 **주파수가 완전히 동일하고 chirality 만 반대**다
(m_polar 가 Pca2₁+시간역전 궤도 안에서 k 를 제자리로 보내고 χ 만 뒤집는다).
따라서 편극 비트는 일반 포논 분광에는 안 보이고 **chirality 감응 측정에만** 보인다.
이는 parent/mirror FORCE_SETS 검증 기준이기도 하다 — 둘의 주파수는 모든 k 에서
같아야 하며, 레포가 관측한 0.14 THz 불일치는 전적으로 DFT 설정 차이다.
"""

import argparse
import itertools
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph                       # noqa: E402

BAND = 17
LOW = BAND - 1


def chirality(ph, kc, R=0.004, nt=24, nph=48):
    th = np.linspace(0, np.pi, nt)
    phi = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    K = np.array([[kc + R * np.array([np.sin(t) * np.cos(p),
                                      np.sin(t) * np.sin(p),
                                      np.cos(t)]) for p in phi] for t in th]).reshape(-1, 3)
    ph.run_qpoints(K.tolist(), with_eigenvectors=True)
    ev = ph.qpoints.eigenvectors
    u = np.array([ev[n][:, LOW] for n in range(len(K))]).reshape(nt, nph, -1)
    tot = 0.0
    for i in range(nt - 1):
        for j in range(nph):
            j2 = (j + 1) % nph
            a, b, c, d = u[i, j], u[i, j2], u[i + 1, j2], u[i + 1, j]
            z = np.vdot(a, b) * np.vdot(b, c) * np.vdot(c, d) * np.vdot(d, a)
            if abs(z) > 0:
                tot += np.angle(z)
    return tot / (2 * np.pi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='source/pristine_mirror')
    ap.add_argument('--seed', type=float, nargs=3, default=[0.097523, 0.0, 0.161028],
                    help='Weyl 하나의 좌표. 나머지는 부호 궤도로 만든다.')
    ap.add_argument('--polar', type=int, default=1)
    ap.add_argument('--dim', type=int, nargs=3, default=[2, 2, 2])
    a = ap.parse_args()

    ph = get_ph(a.dir, tuple(a.dim))
    k0 = np.array(a.seed)
    free = [i for i in range(3) if i != a.polar]

    print('Weyl 궤도 chirality  (%s, band %d-%d)' % (a.dir, BAND, BAND + 1))
    orb = []
    for signs in itertools.product((+1, -1), repeat=2):
        k = k0.copy()
        for i, s in zip(free, signs):
            k[i] = s * k0[i]
        ph.run_qpoints([list(k)])
        f = ph.qpoints.frequencies[0]
        c = chirality(ph, k)
        orb.append((k.copy(), c))
        print('   k=(%+.6f, %+.6f, %+.6f)  gap=%.2e  chi=%+.3f'
              % (k[0], k[1], k[2], f[BAND] - f[LOW], c))

    print('\n표면/계면 BZ 투영별 알짜 chirality  (0 이 아니면 arc 존재)')
    ax = 'abc'
    for nrm in range(3):
        tag = ' (편극축)' if nrm == a.polar else ''
        proj = {}
        for k, c in orb:
            key = tuple(round(float(k[i]), 6) for i in range(3) if i != nrm)
            proj[key] = proj.get(key, 0.0) + c
        arc = any(abs(v) > 0.5 for v in proj.values())
        print('  법선 ∥ %s%s -> 투영면 (%s)   %s'
              % (ax[nrm], tag, ','.join('k%d' % (i + 1) for i in range(3) if i != nrm),
                 'ARC 있음' if arc else 'arc 없음 (알짜 0)'))
        for key, v in sorted(proj.items()):
            print('       투영점 (%+.6f, %+.6f) : 알짜 chi = %+.2f' % (key[0], key[1], v))


if __name__ == '__main__':
    main()
