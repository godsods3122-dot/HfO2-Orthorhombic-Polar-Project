#!/usr/bin/env python3
"""
polar_plane_track.py — 편극평면(k_P = 0) 안의 노드를 strain 계열에 걸쳐 추적한다.

**축 규약 문제를 우회한다.** 구조마다 어느 지표가 편극축인지가 다르므로
(CLAUDE.md §2.1), 이 스크립트는 spglib 에서 편극축 지표를 직접 뽑고, 남은 두
거울축을 격자상수로 **ML(긴 거울축) / MS(짧은 거울축)** 로 이름 붙인다.
좌표를 전부 이 물리 라벨 `(k_ML, k_MS, k_P)` 로 출력하므로 구조끼리 바로 비교된다.

탐색은 CLAUDE.md 의 경고를 지킨다 — gap 순위로 자르지 않고 격자 국소최소를
**전부** 3차원 자유 미세화(구속 없음)하고, 판정은 도달점이 일반위치에 머무는지로만
한다. 3차원 자유 미세화라서 노드가 편극평면을 떠나면 그것도 보인다.

사용법
------
    python3 scripts/polar_plane_track.py <dir> <band> [N] [thresh] [--seed kML kMS kP]

    python3 scripts/polar_plane_track.py source/m1_mirror 18 161 0.03
    python3 scripts/polar_plane_track.py source/m1_mirror 18 0 0 --seed 0.3014 0.0217 0.0

N=0 이면 격자 스캔을 건너뛰고 --seed 만 미세화한다 (계열 연속 추적용).
--band 18 이면 18번과 19번 밴드 사이 gap 을 본다.
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph, refine, symmetry_info


def axis_labels(d):
    """편극축 지표와 ML/MS 지표를 0-based 로 돌려준다."""
    from phonopy.interface.calculator import read_crystal_structure
    name, mirrors, polar, _, _ = symmetry_info(d)
    if len(polar) != 1:
        sys.exit("편극축 지표가 %s 로 하나가 아니다. 수동 확인 필요." % polar)
    p = polar[0]
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    L = np.linalg.norm(u.cell, axis=1)
    ml, ms = sorted(mirrors, key=lambda i: -L[i])
    return name, p, ml, ms, L


def main():
    argv = sys.argv[1:]
    seed = None
    if '--seed' in argv:
        i = argv.index('--seed')
        seed = [float(x) for x in argv[i+1:i+4]]
        argv = argv[:i] + argv[i+4:]
    if len(argv) < 2:
        sys.exit(__doc__)
    d = argv[0]
    b = int(argv[1])
    N = int(argv[2]) if len(argv) > 2 else 161
    thr = float(argv[3]) if len(argv) > 3 else 0.03

    name, ip, iml, ims, L = axis_labels(d)
    print("%s  %s" % (d, name))
    print("  편극축 P = k%d (%.5f Å)   ML = k%d (%.5f Å)   MS = k%d (%.5f Å)"
          % (ip + 1, L[ip], iml + 1, L[iml], ims + 1, L[ims]))

    def native(kml, kms, kp):
        q = [0.0, 0.0, 0.0]
        q[iml] = kml; q[ims] = kms; q[ip] = kp
        return q

    def phys(q):
        return (q[iml], q[ims], q[ip])

    ph = get_ph(d, (2, 2, 2))
    B = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi
    bl = np.linalg.norm(B, axis=1)

    def report(tag, q, g):
        kml, kms, kp = phys(q)
        ph.run_qpoints([list(q)])
        f = ph.qpoints.frequencies[0]
        onmir = (min(abs(kml), abs(abs(kml) - .5)) < 3e-4 or
                 min(abs(kms), abs(abs(kms) - .5)) < 3e-4)
        onnp = min(abs(abs(kp) - .5), 1) < 3e-4
        tagg = "거울면(χ=0)" if onmir else ("nodal plane(χ=0)" if onnp else "일반위치")
        print("  %s (k_ML,k_MS,k_P)=(%.7f,%.7f,%.7f) gap=%.3e  E=%.6f THz = %.3f meV"
              % (tag, kml, kms, kp, g, f[b-1], f[b-1] * 4.135667696))
        print("      MS 거울면까지 %.5f 환산 = %.5f 1/Å   |  편극평면에서 %.5f 환산   [%s]"
              % (abs(kms), abs(kms) * bl[ims], abs(kp), tagg))
        return tagg

    if seed is not None:
        x, g = refine(ph, native(*seed), b)
        print("\n--- seed 연속 추적 (구속 없는 3D refine) ---")
        report("seed->", list(x), g)

    if N and N > 2:
        us = np.linspace(0.0, 0.5, N)
        G = np.zeros((N, N))
        for i, a in enumerate(us):
            ph.run_qpoints([native(a, c, 0.0) for c in us])
            f = ph.qpoints.frequencies
            G[i] = f[:, b] - f[:, b-1]
            if i % 40 == 0:
                print("  k_ML %d/%d" % (i, N), flush=True)
        print("\n편극평면 최소 gap = %.4e  (band %d-%d)" % (G.min(), b, b+1))
        cands = []
        for i in range(1, N-1):
            for j in range(1, N-1):
                v = G[i, j]
                if v < thr and v == G[i-1:i+2, j-1:j+2].min():
                    cands.append((v, us[i], us[j]))
        cands.sort()
        print("gap<%.3f 격자 국소최소 %d개 -> 전부 3D 자유 미세화" % (thr, len(cands)))
        seen, gen = [], 0
        for v, a, c in cands:
            x, g = refine(ph, native(a, c, 0.0), b)
            key = tuple(np.round(np.abs(phys(list(x))), 4))
            if key in seen:
                continue
            seen.append(key)
            t = report("", list(x), g)
            if t == "일반위치":
                gen += 1
        print("\n  서로 다른 도달점 %d개, 그중 일반위치 %d개" % (len(seen), gen))


if __name__ == '__main__':
    main()
