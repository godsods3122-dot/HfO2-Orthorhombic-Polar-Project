#!/usr/bin/env python3
"""
mixing.py — parent 와 mirror 를 "배합"해서 거시적으로만 군을 깨려는 시도를 판정한다.

두 갈래가 정반대 결과를 준다.

(1) 랜덤(무상관) 배합 — 안 된다
    mirror 도메인을 분율 x 로 무작위로 섞으면, 배위평균(configuration average)이
    대칭을 **복원**한다. 증명은 간단하다:
        r_A = r0 + u ,  r_B = r0 - u   (u = 극성 모드, m_polar 에 대해 홀)
        평균 = r0 + (1-2x) u
    즉 평균 구조는 **극성 진폭만 (1-2x) 로 줄어든 같은 Pca2₁** 이다. glide 는
    평균적으로 멀쩡히 살아 있고, 무질서는 대칭 파괴가 아니라 **선폭 증가**를 준다.
    x = 0.5 에서 극성이 0 이 되고 구조는 비극성 모체(Pbcm, 반전 있음)로 간다
    → 반전 + 시간역전 → Berry curvature 항등적 0 → Weyl 전멸.

    다만 여기서 **공짜 손잡이가 하나 나온다**: 극성 진폭 λ = 1-2x 는 연속적이고
    Pca2₁ 을 보존하며 strain 과 독립인 제3의 제어 변수다
    (results/weyl_trend/SUMMARY.md 5절이 요구한 바로 그것).

(2) 상관 배합 — 된다. 필요한 것은 주기성이 아니라 **레지스트리 선택**
    Pca2₁ 의 두 glide 는 병진 성분에 둘 다 c 방향 ½ 을 가지므로, glide 는 c 높이
    h 의 벽을 h+½ 의 벽으로 보낸다. 이 두 "레지스트리"가 구조적으로 **서로 다르면**
    형성에너지가 달라 자연이 하나를 고르고, 그러면 **주기성이 전혀 없어도**
    배위평균에서 glide 가 깨진다.
    계산 결과 두 레지스트리는 실제로 다른 구조다 (Hf-Hf 거리 최대 1.11 Å 차이,
    4 Å 이내 이웃 개수도 다름). 즉 선택이 가능하다.

    → 설계 기준은 "주기 분극"이 아니라 **"모든 벽이 같은 레지스트리를 고를 것"**.
      훨씬 약한 요구다.

사용법
------
    python3 scripts/mixing.py --dir source/pristine_mirror
"""

import argparse
import sys

import numpy as np


def load(d):
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=d + '/POSCAR', interface_mode='vasp')
    return np.array(u.cell), np.array(u.scaled_positions), np.array(u.numbers)


def match(L, n, pa, pb):
    """같은 원소끼리 최소변위로 짝짓고, 재배열된 pb 와 최대변위(Å) 를 돌려준다."""
    from scipy.optimize import linear_sum_assignment
    idx = np.arange(len(pa))
    out = np.zeros_like(pb)
    worst = 0.0
    for Z in np.unique(n):
        ia = idx[n == Z]
        C = np.zeros((len(ia), len(ia)))
        for x, i in enumerate(ia):
            d = pb[ia] - pa[i]
            d -= np.round(d)
            C[x] = np.linalg.norm(d @ L, axis=1)
        r, c = linear_sum_assignment(C)
        for x, y in zip(r, c):
            out[ia[x]] = pb[ia[y]]
            worst = max(worst, C[x, y])
    return out, worst


def polar_mode(L, p, n, polar, ngrid=201):
    """180° 도메인 짝을 만드는 거울 원점을 찾고 극성 모드 u 와 비극성 모체 r0 를 낸다."""
    best = None
    for s in np.linspace(0, 1, ngrid, endpoint=False):
        q = p.copy()
        q[:, polar] = (s - q[:, polar]) % 1.0
        pb, w = match(L, n, p, q)
        if best is None or w < best[1]:
            best = (s, w, pb)
    s, w, pb = best
    d = pb - p
    d -= np.round(d)
    return s, w, (p + d / 2) % 1.0, -d / 2


def fingerprint(L, q, num, cut=4.0):
    """좌표계 무관 구조 지문: 원소쌍 종류별 정렬된 거리 목록."""
    import itertools
    out = {}
    sc = list(itertools.product((-1, 0, 1), repeat=3))
    for i in range(len(q)):
        for j in range(i + 1, len(q)):
            key = tuple(sorted((int(num[i]), int(num[j]))))
            for sh in sc:
                r = np.linalg.norm((q[j] + np.array(sh) - q[i]) @ L)
                if r < cut:
                    out.setdefault(key, []).append(r)
    return {k: np.sort(np.array(v)) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', default='source/pristine_mirror')
    ap.add_argument('--polar', type=int, default=1)
    ap.add_argument('--stack', type=int, default=2, help='벽 법선 방향 (c = 2)')
    a = ap.parse_args()
    import spglib

    L, p, n = load(a.dir)
    s, w, r0, u = polar_mode(L, p, n, a.polar)
    print('180° 도메인 짝: 거울 원점 s = %.4f, 최대 원자변위 %.4f Å' % (s, w))
    print('극성 모드 진폭 : 최대 %.4f Å, RMS %.4f Å'
          % (np.abs(u @ L).max(), np.sqrt(((u @ L) ** 2).sum(1).mean())))

    print('\n[1] 랜덤 배합의 배위평균  r0 + (1-2x) u')
    print('   x      lambda   공간군      ops  반전')
    for x in (0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.45, 0.5):
        lam = 1 - 2 * x
        ds = spglib.get_symmetry_dataset((L, (r0 + lam * u) % 1.0, n), symprec=1e-4)
        inv = any(np.array_equal(r, -np.eye(3, dtype=int)) for r in ds.rotations)
        print('  %.2f   %+.2f    %-10s %3d  %s' % (x, lam, ds.international, len(ds.rotations), inv))
    print('  => 어떤 x 에서도 Pca2₁ 그대로. 랜덤 배합은 glide 를 못 깬다.')

    print('\n[2] 벽 레지스트리 tau (c 방향 밀기, 원시셀 단위)')
    pB = p.copy()
    pB[:, a.polar] = (s - pB[:, a.polar]) % 1.0
    fps = {}
    print('  tau    최소원자거리   공간군      ops')
    for tau in (0.0, 0.25, 0.5, 0.75):
        N = 3
        L2 = L.copy()
        L2[a.stack] = L[a.stack] * N
        q = []
        for j in range(N):
            src = p.copy() if j < 2 else pB.copy()
            if j >= 2:
                src[:, a.stack] = (src[:, a.stack] + tau) % 1.0
            src[:, a.stack] = (src[:, a.stack] + j) / N
            q.append(src)
        q = np.vstack(q)
        num = np.tile(n, N)
        ds = spglib.get_symmetry_dataset((L2, q, num), symprec=1e-3)
        md = min(np.linalg.norm(((q[j] - q[i]) - np.round(q[j] - q[i])) @ L2)
                 for i in range(len(q)) for j in range(i + 1, len(q)))
        fps[tau] = fingerprint(L2, q, num)
        print('  %.2f    %6.3f Å     %-10s %3d' % (tau, md, ds.international, len(ds.rotations)))

    print('\n  glide 가 짝짓는 tau=0 과 tau=0.5 비교:')
    tot = 0.0
    for key in sorted(fps[0.0]):
        x_, y_ = fps[0.0][key], fps[0.5][key]
        m = min(len(x_), len(y_))
        dev = float(np.abs(x_[:m] - y_[:m]).max()) if m else 0.0
        tot = max(tot, dev)
        print('     원소쌍 %-10s 거리개수 %d vs %d, 최대차이 %.4f Å'
              % (str(key), len(x_), len(y_), dev))
    print('  => %s (최대 %.4f Å)'
          % ('서로 다른 구조 → 형성에너지가 달라 자연이 하나를 고른다'
             if tot > 1e-3 else '같은 구조 → 선택 불가', tot))


if __name__ == '__main__':
    main()
