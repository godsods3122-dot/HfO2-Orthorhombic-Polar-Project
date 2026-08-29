#!/usr/bin/env python3
"""
order_break.py — Pca2_1 의 glide 두 개만 부분적으로 깨서(→ P2_1) 거울면 nodal line을
고대칭 경로 위의 진짜 Weyl point로 바꾸는 경로를 검증한다.

배경
----
Pca2_1(mm2)에서 chirality=0 을 강제하는 것은 **glide 거울 두 개**뿐이다.
2_1 나사축(= 편극축, 강유전성의 근원)은 k_polar=0.5 nodal plane 하나만 강제하고
고대칭 경로 위에는 아무 제약도 걸지 않는다. 따라서

    mm2 (Pca2_1)  ->  2 (P2_1)      [glide 2개 제거, 2_1 유지]

는 **강유전성과 편극축을 그대로 두면서** 표준 고대칭 경로를 여는 최소 대칭 축소다.
shear/전단은 전혀 필요 없다 — 조성·미세구조 규칙배열로 도달한다.

핵심 관찰
---------
거울면(k1=0) 안에는 band17-18 nodal line 이 있고, 그것이 Gamma-Z 를 k3~0.2185 에서
가로지른다. Pca2_1 에서는 거울이 chi=0 을 강제하므로 위상적으로 무의미하다.
glide 를 깨면 이 nodal line 은 gap 이 열리지만, k_polar=0 평면 위의 한 점만은
C2*T (반유니타리, 각 k 를 고정 -> 국소적으로 실수 해밀토니안) 로 보호되어 살아남고,
거울이 없어졌으므로 이제 **chi=+-1 인 진짜 Weyl** 이 된다.

남는 문제는 그 점이 glide-breaking 세기에 비례해 k1=0 에서 벗어난다는 것인데,
**서로 독립인 규칙배열 채널 2개**(양이온 / 산소)를 비율로 섞으면 k1=0 을 정확히
유지한 채 glide 만 깰 수 있다. 조건이 1개(k1=0), 자유도가 1개(채널 비율)이므로
항상 풀린다.

이 스크립트가 하는 일
--------------------
  --mode enum     4a 궤도 안에서 어떤 2원자 치환이 P2_1 을 주는지 spglib 으로 열거
  --mode response 각 규칙배열 채널의 dk1/d(delta) 측정
  --mode solve    두 채널을 섞어 Weyl 을 Gamma-Z 위에 정확히 올리고 chirality 계산

주의 — 이 스크립트의 규칙배열은 **질량만** 바꾼다. 따라서
  * 동위원소 규칙배열(16O/18O 등)에 대해서는 Born-Oppenheimer 상 **정확**하다.
  * Hf/Zr 화학 치환에 대해서는 힘상수·Born charge 변화를 뺀 하한 모형이다
    (실제 효과는 더 크다). 메커니즘 자체는 대칭 논증이라 이 근사와 무관하다.
"""

import argparse
import itertools
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph                      # noqa: E402

BAND = 17                      # band17-18 gap (0-based 16,17)
LOW = BAND - 1


# ----------------------------------------------------------------- 대칭 열거
def enumerate_orderings(d, symprec=1e-4):
    """4a 궤도를 2+2 로 쪼개는 모든 방법의 공간군을 돌려준다."""
    import spglib
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    ds = spglib.get_symmetry_dataset((u.cell, u.scaled_positions, u.numbers), symprec=symprec)
    orbits = {}
    for i, e in enumerate(ds.equivalent_atoms):
        orbits.setdefault(int(e), []).append(i)
    out = []
    for key, orb in orbits.items():
        for pair in itertools.combinations(orb, 2):
            n = u.numbers.copy()
            n[list(pair)] = 200 - key            # 궤도마다 다른 더미 원소
            sg = spglib.get_symmetry_dataset((u.cell, u.scaled_positions, n),
                                             symprec=symprec).international
            out.append((key, orb, pair, sg))
    return ds.international, orbits, out


# --------------------------------------------------------------- 질량 규칙배열
def apply_channels(ph, m0, channels):
    """channels: [(heavy_pair, light_pair, delta), ...]  — 평균질량 보존."""
    m = np.array(m0, float)
    for hp, lp, dl in channels:
        m[list(hp)] += dl
        m[list(lp)] -= dl
    ph.masses = list(m)


def gaps(ph, qs, with_freq=False):
    ph.run_qpoints([list(q) for q in qs])
    f = ph.qpoints.frequencies
    g = f[:, BAND] - f[:, LOW]
    return (g, f[:, LOW]) if with_freq else g


def refine_in_polar_plane(ph, seed, polar=1):
    """편극축 지표를 0 으로 고정한 평면(C2*T 보호) 안에서 노드를 정밀화."""
    from scipy.optimize import minimize
    ax = [i for i in range(3) if i != polar]

    def f(x):
        q = [0.0, 0.0, 0.0]
        q[ax[0]], q[ax[1]] = x
        return float(gaps(ph, [q])[0])

    r = minimize(f, seed, method='Nelder-Mead',
                 options=dict(xatol=1e-12, fatol=1e-18, maxiter=9000, maxfev=9000))
    return r.x, r.fun


def chirality(ph, kc, R=0.004, nt=22, nph=44):
    """band LOW 의 Berry flux 를 구면에서 Fukui-Hatsugai-Suzuki 로 적분.

    Simphony 의 WeylChirality_calc 와 전체 부호 규약이 반대다(크기·상대부호는 일치).
    """
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
            a, b, c, dd = u[i, j], u[i, j2], u[i + 1, j2], u[i + 1, j]
            z = np.vdot(a, b) * np.vdot(b, c) * np.vdot(c, dd) * np.vdot(dd, a)
            if abs(z) > 0:
                tot += np.angle(z)
    return tot / (2 * np.pi)


# --------------------------------------------------------------------- main
CATION = ((0, 3), (1, 2))          # pristine_mirror 의 C2 짝 (op1: 0<->3, 1<->2)
OXYGEN = ((4, 11), (5, 10))        # O 궤도 1 의 C2 짝
SEED = [0.0, 0.218545]             # Gamma-Z 위 nodal line 교차점 (k1, k3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='source/pristine_mirror')
    ap.add_argument('--mode', default='solve',
                    choices=['enum', 'response', 'solve'])
    ap.add_argument('--dim', type=int, nargs=3, default=[2, 2, 2])
    ap.add_argument('--polar', type=int, default=1, help='편극축 지표 (0-based)')
    ap.add_argument('--dhf', type=float, default=4.0, help='양이온 채널 세기 (amu)')
    a = ap.parse_args()

    if a.mode == 'enum':
        sg, orbits, rows = enumerate_orderings(a.dir)
        print('parent space group: %s' % sg)
        print('4a orbits: %s' % orbits)
        print('\n2+2 치환 -> 공간군  (P2_1 인 것만이 glide 2개를 지우고 2_1 을 남긴다)')
        for key, orb, pair, s in rows:
            mark = '  <== P2_1' if s == 'P2_1' else ''
            print('   orbit %2d  치환 %-10s -> %-8s%s' % (key, str(list(pair)), s, mark))
        return

    ph = get_ph(a.dir, tuple(a.dim))
    m0 = np.array(ph.masses, float)

    if a.mode == 'response':
        print('규칙배열 채널별 dk1/d(delta)   (Gamma-Z 위 노드 %s 기준)' % SEED)
        for lab, ch in [('cation', CATION), ('oxygen', OXYGEN)]:
            for dl in (+1.0, -1.0):
                apply_channels(ph, m0, [(ch[0], ch[1], dl)])
                x, g = refine_in_polar_plane(ph, SEED, a.polar)
                print('   %-7s delta=%+.1f amu : k=(%+.7f, %+.7f) gap=%.2e  dk1=%+.7f'
                      % (lab, dl, x[0], x[1], g, x[0] - SEED[0]))
        return

    # ---- solve : 두 채널을 섞어 k1 = 0 을 유지한 채 glide 를 깬다
    from scipy.optimize import brentq
    apply_channels(ph, m0, [(CATION[0], CATION[1], 1.0)])
    aC = refine_in_polar_plane(ph, SEED, a.polar)[0][0]
    apply_channels(ph, m0, [(OXYGEN[0], OXYGEN[1], 1.0)])
    aO = refine_in_polar_plane(ph, SEED, a.polar)[0][0]
    print('dk1/d(cation) = %+.7f /amu' % aC)
    print('dk1/d(oxygen) = %+.7f /amu' % aO)
    guess = -aC / aO * a.dhf
    print('상쇄 비율 oxygen/cation = %+.6f  ->  dhf=%.4f 이면 dO=%.6f amu'
          % (-aC / aO, a.dhf, guess))

    def k1of(t):
        apply_channels(ph, m0, [(CATION[0], CATION[1], a.dhf), (OXYGEN[0], OXYGEN[1], t)])
        return refine_in_polar_plane(ph, SEED, a.polar)[0][0]

    dO = brentq(k1of, guess - 0.4 * abs(guess) - 0.05, guess + 0.4 * abs(guess) + 0.05, xtol=1e-10)
    apply_channels(ph, m0, [(CATION[0], CATION[1], a.dhf), (OXYGEN[0], OXYGEN[1], dO)])
    x, g = refine_in_polar_plane(ph, SEED, a.polar)
    k = [0.0, 0.0, 0.0]
    ax = [i for i in range(3) if i != a.polar]
    k[ax[0]], k[ax[1]] = x
    gg, fr = gaps(ph, [k], with_freq=True)
    print('\n=== dhf=%.4f amu, dO=%.6f amu ===' % (a.dhf, dO))
    print('  Weyl k = (%+.9f, %+.9f, %+.9f)' % tuple(k))
    print('  gap    = %.3e THz     mode freq = %.5f THz (%.2f cm^-1, %.2f meV)'
          % (g, fr[0], fr[0] * 33.35641, fr[0] * 4.135667))
    print('  chirality = %+.4f' % chirality(ph, np.array(k)))


if __name__ == '__main__':
    main()
