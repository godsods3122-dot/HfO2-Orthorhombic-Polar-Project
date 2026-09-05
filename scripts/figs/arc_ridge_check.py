#!/usr/bin/env python3
"""마스크 인공물 가려내기: 간극 내 최대점이 마스크 경계에서 몇 칸 떨어져 있나.

간극만 남긴 아크 그림에서 끝점을 읽으면 안 되는 이유가 있다 — 노드가 사영 벌크
연속체의 경계에 정확히 앉아 있어서(SUMMARY 6절), 경계에 붙어 달리는 신호는 마스크에
잘려 "노드에서 끝나는" 것처럼 보인다.

가려내는 기준: k_a>0 반쪽에서 행별 간극 내 최대점을 찾아, 그것이 마스크 경계
(그 행의 첫 간극 격자점)에서 몇 칸 떨어져 있는지 본다.
  거리 0  -> 벌크 영역 구조가 잘린 단면. 표면 상태로 볼 근거 없음.
  거리 >0 -> 진짜 간극 내 능선.

사용: python3 arc_ridge_check.py ARC.npz GAPMASK.npy A B
      (A, B = 노드 사영 좌표 |k_a|, |k_b|)
"""
import sys
import numpy as np

npz, maskf = sys.argv[1], sys.argv[2]
A, B = float(sys.argv[3]), float(sys.argv[4])
d = np.load(npz); ks = d['ks']
gap = np.load(maskf)                       # True = gapped, [ka, kb]
Z = np.abs(d['raw_d'])
dk = ks[1] - ks[0]

ing, outg = Z[gap], Z[~gap]
print('|raw_d|   간극 영역: max %8.3f  mean %7.3f  (n=%d)' % (ing.max(), ing.mean(), ing.size))
print('|raw_d| 연속체 영역: max %8.3f  mean %7.3f  (n=%d)' % (outg.max(), outg.mean(), outg.size))
print('노드 사영 (%.4f, %.4f),  격자 dk=%.4f' % (A, B, dk))
print()
print('   k_b    경계 k_a   능선 k_a    세기   경계로부터(칸)')
for j, kb in enumerate(ks):
    if kb < -1e-9 or kb > 2.2 * B:
        continue
    g = gap[:, j] & (ks > 0)
    if not g.any():
        print('  %+6.3f   (이 행에 간극 없음)' % kb)
        continue
    idx = np.where(g)[0]
    edge = idx[0]
    ig = idx[np.argmax(Z[idx, j])]
    tag = '  <- 노드행' if abs(kb - B) < dk / 2 + 1e-9 else ''
    print('  %+6.3f   %+7.3f   %+7.3f   %6.3f   %d%s'
          % (kb, ks[edge], ks[ig], Z[ig, j], ig - edge, tag))
