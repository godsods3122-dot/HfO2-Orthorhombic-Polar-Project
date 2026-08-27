#!/usr/bin/env python3
"""
Nodes.dat에서 nodal line과 고립점(Weyl point 후보)을 분리한다.

원본 대비 바뀐 점
-----------------
1) query_pairs를 파이썬 set이 아니라 ndarray로 받는다.
   8만 점/dist_tol=0.03이면 쌍이 약 2천만 개인데, set은 튜플 객체라
   2GB 이상을 먹고 죽는다. ndarray는 같은 내용이 0.3GB다.

2) 그 전에 격자 스냅으로 점을 솎아낸다.
   FindNodes는 같은 노드를 여러 시작점에서 중복 수렴시키고 노달라인을
   촘촘히 훑기 때문에 대부분이 중복이다. dist_tol보다 작은 격자로 스냅하면
   연결성은 보존하면서 점 수와 쌍 수가 크게 줄어든다.

3) graph + graph.T 를 하지 않는다.
   connected_components(directed=False)가 이미 비대칭 행렬을 무향으로
   해석하므로 대칭화는 메모리만 두 배로 쓴다.

4) 분류 기준을 "점 개수"에서 "공간적 크기"로 바꿨다.
   개수 기준은 오분류가 난다. 진짜 Weyl point도 여러 시작점에서 중복
   수렴하면 점이 수십 개가 되고, 짧은 노달라인 조각은 점이 적을 수 있다.
   고립점은 모든 방향으로 퍼짐이 0에 가깝고 노달라인은 한 방향으로 길다는
   성질을 PCA로 직접 본다.

5) k공간 주기성과 이방성을 옵션으로 처리한다.
   기약좌표에서 0.5와 -0.5는 같은 점인데 KD-Tree는 모른다. periodic=True면
   [0,1)로 감싸고 boxsize=1로 처리한다. cartesian=True면 역격자로 변환해
   실제 거리로 재는데, 축 길이가 다른 사방정계에서 의미가 있다.

사용법
------
    python extract_weyl_fast.py Nodes.dat
    python extract_weyl_fast.py Nodes.dat --dist-tol 0.02 --point-size 0.004
    python extract_weyl_fast.py Nodes.dat --poscar POSCAR --cartesian
"""

import argparse
import sys
import time

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def read_reciprocal(poscar):
    """POSCAR에서 역격자(2*pi 없음, 행이 b1,b2,b3)를 읽는다."""
    with open(poscar) as f:
        lines = f.readlines()
    scale = float(lines[1].split()[0])
    cell = np.array([[float(x) for x in lines[2 + i].split()[:3]] for i in range(3)]) * scale
    return np.linalg.inv(cell).T


def cluster_nodes(k, dist_tol, snap=None, periodic=False):
    """
    격자 스냅으로 대표점을 뽑고, 대표점끼리 연결해서 클러스터를 만든 뒤
    원본 점 전체에 라벨을 되돌려준다.

    반환: labels(원본 길이), n_clusters
    """
    n = len(k)
    if snap is None:
        snap = dist_tol / 3.0

    # --- 격자 스냅: 같은 칸에 떨어진 점들을 대표점 하나로 묶는다 ---
    kk = np.mod(k, 1.0) if periodic else k
    cell_id = np.floor(kk / snap).astype(np.int64)
    _, rep_index, inverse = np.unique(cell_id, axis=0, return_index=True, return_inverse=True)
    inverse = inverse.ravel()
    rep = kk[rep_index]
    m = len(rep)

    # --- 대표점끼리만 쌍을 구한다 ---
    if periodic:
        tree = cKDTree(np.mod(rep, 1.0), boxsize=1.0)
    else:
        tree = cKDTree(rep)
    pairs = tree.query_pairs(r=dist_tol, output_type='ndarray')

    if len(pairs):
        g = coo_matrix((np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])),
                       shape=(m, m))
    else:
        g = coo_matrix((m, m))

    # directed=False면 대칭화 없이도 무향으로 해석한다
    n_rep_clusters, rep_labels = connected_components(g, directed=False)

    # 대표점 라벨을 원본 점으로 전파
    labels = rep_labels[inverse]
    return labels, n_rep_clusters, m, len(pairs)


def classify(k, labels, n_clusters, point_size):
    """클러스터마다 공간적 퍼짐을 재서 고립점/노달라인으로 나눈다."""
    order = np.argsort(labels, kind='stable')
    sorted_labels = labels[order]
    bounds = np.searchsorted(sorted_labels, np.arange(n_clusters + 1))

    isolated, extended = [], []
    for c in range(n_clusters):
        idx = order[bounds[c]:bounds[c + 1]]
        pts = k[idx]
        if len(pts) == 1:
            extent = 0.0
        else:
            centered = pts - pts.mean(axis=0)
            # 가장 긴 주축 방향의 퍼짐 (전체 길이 기준)
            extent = float(np.ptp(centered @ np.linalg.svd(centered, full_matrices=False)[2][0]))
        (isolated if extent <= point_size else extended).append((c, idx, extent))
    return isolated, extended


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('filename', nargs='?', default='Nodes.dat')
    ap.add_argument('--dist-tol', type=float, default=0.03,
                    help='이 거리 안이면 같은 덩어리로 본다 (기약좌표, 기본 0.03)')
    ap.add_argument('--point-size', type=float, default=None,
                    help='이 크기 이하로 뭉쳐 있으면 고립점 (기본 dist-tol/3)')
    ap.add_argument('--snap', type=float, default=None,
                    help='솎아내기 격자 크기 (기본 dist-tol/3)')
    ap.add_argument('--periodic', action='store_true', help='k공간 주기성 적용')
    ap.add_argument('--poscar', default=None)
    ap.add_argument('--cartesian', action='store_true',
                    help='기약좌표 대신 직교 k로 거리를 잰다 (--poscar 필요)')
    ap.add_argument('--max-print', type=int, default=50)
    args = ap.parse_args()

    point_size = args.point_size if args.point_size is not None else args.dist_tol / 3.0

    t0 = time.time()
    try:
        data = np.loadtxt(args.filename, comments='#')
    except OSError:
        sys.exit(f"오류: {args.filename} 을 찾을 수 없습니다.")
    if data.ndim == 1:
        data = data[None, :]
    print(f"[1/4] 로드 {len(data):,}행  ({time.time()-t0:.1f}s)")

    k = data[:, 5:8].copy()          # k1,k2,k3 (기약좌표)
    metric_k = k
    if args.cartesian:
        if not args.poscar:
            sys.exit("--cartesian 을 쓰려면 --poscar 도 주세요.")
        B = read_reciprocal(args.poscar)
        metric_k = k @ B             # 직교 k (1/Angstrom, 2*pi 없음)
        if args.periodic:
            print("      주의: --cartesian 과 --periodic 은 함께 쓰지 않는 것이 안전합니다.")

    t0 = time.time()
    labels, n_clusters, n_rep, n_pairs = cluster_nodes(
        metric_k, args.dist_tol, snap=args.snap, periodic=args.periodic)
    print(f"[2/4] 격자 스냅 {len(k):,} -> 대표점 {n_rep:,}개, 쌍 {n_pairs:,}개  ({time.time()-t0:.1f}s)")
    print(f"[3/4] 클러스터 {n_clusters:,}개")

    isolated, extended = classify(metric_k, labels, n_clusters, point_size)
    print(f"[4/4] 고립점 {len(isolated):,}개 / 늘어난 덩어리(노달라인) {len(extended):,}개\n")

    if extended:
        ext = sorted(extended, key=lambda x: -x[2])[:5]
        print("가장 긴 덩어리 5개 (노달라인으로 판정):")
        for c, idx, e in ext:
            print(f"   점 {len(idx):>6,}개  최대 퍼짐 {e:.4f}  E={data[idx,4].mean():.4f} THz")
        print()

    if not isolated:
        print("고립점이 없습니다. --dist-tol 을 줄이거나 --point-size 를 키워보세요.")
        return

    print(f"고립점 {len(isolated)}개 (클러스터마다 gap이 가장 작은 점 하나씩)")
    print("-" * 92)
    print(f"{'k1':>12} {'k2':>12} {'k3':>12} {'gap':>12} {'E(THz)':>11} {'점수':>6} {'퍼짐':>9}")
    print("-" * 92)
    rows = []
    for c, idx, e in isolated:
        best = idx[np.argmin(data[idx, 3])]
        rows.append((data[best], len(idx), e))
    rows.sort(key=lambda r: r[0][3])
    for pt, cnt, e in rows[:args.max_print]:
        print(f"{pt[5]:12.7f} {pt[6]:12.7f} {pt[7]:12.7f} {pt[3]:12.3e} {pt[4]:11.5f} {cnt:6d} {e:9.5f}")
    if len(rows) > args.max_print:
        print(f"... 외 {len(rows)-args.max_print}개")
    print("-" * 92)

    out = np.array([r[0] for r in rows])
    np.savetxt('isolated_nodes.dat', out,
               header='kx ky kz gap E k1 k2 k3   (gap 오름차순)', fmt='%18.10f')
    print("isolated_nodes.dat 에 저장했습니다.")


if __name__ == '__main__':
    main()
