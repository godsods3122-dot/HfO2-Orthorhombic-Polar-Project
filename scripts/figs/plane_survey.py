#!/usr/bin/env python3
"""polar=0 평면 위 band17/18 노드 전수 조사 — pristine / m1 / p1.

node_audit 의 교훈대로: gap 순위로 자르지 않고 **국소최소 후보를 전부** 격자
미세화한 뒤, (a) 거울선에 빨려들지 않고 일반위치에 남는가 (b) chirality 가 여러
반지름에서 정수로 안정한가 로 판정한다.
"""
import sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0])
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph, gap_at, symmetry_info
from axes import roles

BAND = 17
N = 201                      # 평면 격자 (간격 0.0025)
TOLMIR = 0.012               # 거울선에서 이만큼 안쪽만 후보로

def flux(ph, rec, center, r, nth=36, nph=72, nocc=17):
    th = np.linspace(0, np.pi, nth); phi = np.linspace(0, 2*np.pi, nph, endpoint=False)
    T, P = np.meshgrid(th, phi, indexing='ij')
    d = np.stack([np.sin(T)*np.cos(P), np.sin(T)*np.sin(P), np.cos(T)], -1)
    Q = (center + (d*r) @ np.linalg.inv(rec)).reshape(-1, 3)
    ph.run_qpoints([list(q) for q in Q], with_eigenvectors=True)
    V = ph.qpoints.eigenvectors; nb = V.shape[1]
    V = V.reshape(nth, nph, nb, nb)[:, :, :, :nocc]
    L = lambda A, B: (lambda m: m/np.abs(m))(np.linalg.det(
        np.einsum('...ji,...jk->...ik', A.conj(), B)))
    Vp = np.concatenate([V, V[:, :1]], axis=1)
    F = np.angle(L(Vp[:-1,:-1],Vp[:-1,1:]) * L(Vp[:-1,1:],Vp[1:,1:]) /
                 (L(Vp[1:,:-1],Vp[1:,1:]) * L(Vp[:-1,:-1],Vp[1:,:-1])))
    return F.sum()/(2*np.pi)

for src in ('pristine_mirror', 'm1_mirror', 'p1_mirror'):
    D = '/home/user/HfO2-Orthorhombic-Polar-Project/source/' + src
    r2n, pol, mir = roles(D)
    ph = get_ph(D)
    rec = np.linalg.inv(ph.primitive.cell).T * 2*np.pi
    i0, i1 = mir                                     # 거울 지표 (native)
    print('=' * 78); print(src, '  polar native', pol, '  mirror native', mir, flush=True)

    ax = np.linspace(0, 0.5, N)
    A, B = np.meshgrid(ax, ax, indexing='ij')
    Q = np.zeros((N*N, 3)); Q[:, i0] = A.ravel(); Q[:, i1] = B.ravel(); Q[:, pol] = 0.0
    t = time.time(); g = np.empty(N*N)
    for s in range(0, N*N, 20000):
        g[s:s+20000] = gap_at(ph, [list(q) for q in Q[s:s+20000]], BAND)
    g = g.reshape(N, N)
    print('  평면 스캔 %dx%d  %.0fs   전체 최소 %.3e' % (N, N, time.time()-t, g.min()), flush=True)

    # 거울선에서 떨어진 국소최소 후보
    inner = (ax > TOLMIR) & (ax < 0.5 - TOLMIR)
    m = np.where(np.outer(inner, inner), g, np.inf)
    loc = []
    for i in range(1, N-1):
        for j in range(1, N-1):
            if np.isinf(m[i, j]): continue
            if m[i, j] <= m[i-1:i+2, j-1:j+2].min():
                loc.append((m[i, j], ax[i], ax[j]))
    loc.sort()
    print('  거울선 밖 국소최소 %d개 — 전부 미세화한다' % len(loc), flush=True)

    for gv, a0, b0 in loc[:10]:
        ctr = np.array([a0, b0]); half = 0.0025
        for _ in range(7):
            aa = np.linspace(ctr[0]-half, ctr[0]+half, 13)
            bb = np.linspace(ctr[1]-half, ctr[1]+half, 13)
            AA, BB = np.meshgrid(aa, bb, indexing='ij')
            QQ = np.zeros((169, 3)); QQ[:, i0] = AA.ravel(); QQ[:, i1] = BB.ravel(); QQ[:, pol] = 0.0
            gg = gap_at(ph, [list(q) for q in QQ], BAND)
            k = int(np.argmin(gg)); ctr = np.array([QQ[k, i0], QQ[k, i1]]); best = gg[k]; half /= 3
        onmir = min(abs(ctr[0]), abs(ctr[0]-0.5), abs(ctr[1]), abs(ctr[1]-0.5)) < 1e-3
        q = np.zeros(3); q[i0], q[i1] = ctr
        ph.run_qpoints([list(q)]); E = ph.qpoints.frequencies[0][BAND-1]
        tag = '거울선으로 배수' if onmir else '일반위치 유지'
        line = ('    격자 %.3e  ->  native(k%d,k%d)=(%.7f, %.7f)  gap %.3e  E %.6f  [%s]'
                % (gv, i0+1, i1+1, ctr[0], ctr[1], best, E, tag))
        if not onmir and best < 1e-4:
            chis = [flux(ph, rec, q, r) for r in (0.002, 0.004, 0.008)]
            line += '   chi = %s' % ' / '.join('%+.2f' % c for c in chis)
        print(line, flush=True)
