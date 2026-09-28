#!/usr/bin/env python3
"""
box_scan.py — k 공간의 직육면체 박스만 훑어 노드를 찾는다.
"Γ-Z 경로 근처" 처럼 특정 영역을 볼 때 쓴다 (weyl_scan.py 의 wedge 는 쐐기 전체).

CLAUDE.md 의 경고를 그대로 지킨다: **gap 순위로 후보를 자르지 않고 전부 미세화**하고,
"미세화 도달점이 일반위치에 머무는가" 로만 판정한다.

거울 지표는 spglib 에서 뽑지 않고 인자로 준다 (구조마다 축 규약이 다르므로
먼저 `weyl_scan.py --mode sym` 으로 확인할 것).

사용법
------
    python3 scripts/box_scan.py <dir> <band> <k1k2max> <N12> <N3> <thresh> <mirror_idx_0based...>

    # parent_pristine 은 거울 지표가 1-based [1,2] = 0-based 0,1
    python3 scripts/box_scan.py source/parent_pristine 18 0.20 41 51 0.03 0 1

--band 18 이면 18번과 19번 밴드 사이 gap 을 본다.
박스는 k1,k2 in [0, k1k2max], k3 in [0, 0.5] 이다.
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph, refine, on_mirror


def main():
    if len(sys.argv) < 8:
        sys.exit(__doc__)
    d = sys.argv[1]
    b = int(sys.argv[2])
    kmax = float(sys.argv[3])
    n12 = int(sys.argv[4])
    n3 = int(sys.argv[5])
    thr = float(sys.argv[6])
    mir = [int(x) for x in sys.argv[7:]]

    ph = get_ph(d, (2, 2, 2))
    k1s = np.linspace(0, kmax, n12)
    k2s = np.linspace(0, kmax, n12)
    k3s = np.linspace(0, 0.5, n3)
    G = np.zeros((n12, n12, n3), dtype=np.float32)
    for i, a in enumerate(k1s):
        ph.run_qpoints([[a, y, z] for y in k2s for z in k3s])
        f = ph.qpoints.frequencies
        G[i] = (f[:, b] - f[:, b - 1]).reshape(n12, n3)
        if i % 10 == 0:
            print("  k1 %d/%d" % (i, n12), flush=True)

    print("\nband %d-%d, 박스 k1,k2<=%.3f  최소 gap = %.4e" % (b, b + 1, kmax, G.min()))
    cands = []
    for i in range(n12):
        for j in range(n12):
            for k in range(n3):
                v = G[i, j, k]
                if v >= thr:
                    continue
                sl = G[max(0, i-1):i+2, max(0, j-1):j+2, max(0, k-1):k+2]
                if v == sl.min():
                    cands.append((float(v), k1s[i], k2s[j], k3s[k]))
    cands.sort()
    print("gap<%.3f 국소최소 %d개 -> 전부 미세화" % (thr, len(cands)))

    weyl, triv = [], 0
    for v, a, bb, c in cands:
        x, g = refine(ph, [a, bb, c], b)
        if on_mirror(x, mir, tol=3e-4):
            triv += 1
        else:
            weyl.append((g, x))
    print("  거울면으로 배수/거울면 위 : %d개  (chirality 0 강제)" % triv)
    print("  일반 위치에 남은 것       : %d개" % len(weyl))
    for g, x in sorted(weyl, key=lambda t: t[0])[:25]:
        ph.run_qpoints([list(x)])
        f = ph.qpoints.frequencies[0]
        print("   k=(%.7f,%.7f,%.7f) gap=%.3e  E=%.6f THz" % (x[0], x[1], x[2], g, f[b-1]))


if __name__ == '__main__':
    main()
