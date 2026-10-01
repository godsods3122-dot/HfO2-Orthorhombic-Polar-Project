#!/usr/bin/env python3
"""m2.5 Weyl 재검증 — 순수 격자(gradient refine 없음).

gradient refine 은 type-II 콘에서 골을 따라 미끄러져 엉뚱한 데로 수렴한다.
그래서 격자 단계적 미세화만 쓴다.
"""
import sys, time
import numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph, symmetry_info, gap_at

D = sys.argv[1] if len(sys.argv) > 1 else 'source/m2.5_mirror'
BAND = 17
ph = get_ph(D)
name, mir, pol, rots, trans = symmetry_info(D)
print('SG %s   mirror indices %s   polar index %s' % (name, mir, pol), flush=True)
P = pol[0]

def gaps(qs, chunk=20000):
    out = []
    for i in range(0, len(qs), chunk):
        out.append(gap_at(ph, qs[i:i+chunk], BAND))
    return np.concatenate(out)

def freqs(q):
    ph.run_qpoints([list(q)])
    return ph.qpoints.frequencies[0]

# ---- 1. 대칭이 강제하는 평면 확인 -------------------------------------------
print('\n[1] 대칭 강제 평면')
rng = np.random.default_rng(0)
for lab, build in [
    ('polar = 0.5 (nodal plane?)', lambda r: np.insert(r, P, 0.5)),
    ('polar = 0.0',                lambda r: np.insert(r, P, 0.0)),
    ('mirror k%d = 0.0' % (mir[0]+1), lambda r: np.insert(r, mir[0], 0.0)),
    ('mirror k%d = 0.5' % (mir[0]+1), lambda r: np.insert(r, mir[0], 0.5)),
    ('mirror k%d = 0.0' % (mir[1]+1), lambda r: np.insert(r, mir[1], 0.0)),
    ('mirror k%d = 0.5' % (mir[1]+1), lambda r: np.insert(r, mir[1], 0.5)),
]:
    qs = [build(rng.random(2) * 0.5) for _ in range(200)]
    g = gaps(qs)
    print('  %-28s gap  min %.3e  median %.3e  max %.3e' % (lab, g.min(), np.median(g), g.max()))

# ---- 2. 기약 쐐기 전역 격자 --------------------------------------------------
print('\n[2] 기약 쐐기 전역 격자')
for N in (41, 61):
    t = time.time()
    ax = np.linspace(0, 0.5, N)
    Q = np.array(np.meshgrid(ax, ax, ax, indexing='ij')).reshape(3, -1).T
    g = gaps([list(q) for q in Q]).reshape(N, N, N)
    # 대칭 강제 자리 마스킹: 편극축 0.5 평면, 거울면 (지표 0 또는 0.5)
    tol = 1.0 / (N - 1) / 2 * 1e-6 + 1e-9
    K = ax
    on = np.zeros((N, N, N), bool)
    idx = [np.abs(K) < 1e-9, np.abs(K - 0.5) < 1e-9]
    sh = [(-1,1,1),(1,-1,1),(1,1,-1)]
    for i in mir:
        m = (idx[0] | idx[1]).reshape(sh[i])
        on |= np.broadcast_to(m, (N,N,N))
    m = idx[1].reshape(sh[P])
    on |= np.broadcast_to(m, (N,N,N))
    gg = np.where(on, np.inf, g)
    f = np.unravel_index(np.argmin(gg), gg.shape)
    print('  N=%d (%d pts, %.0fs)  일반위치 최소 gap %.4e THz  at (%.4f, %.4f, %.4f)'
          % (N, N**3, time.time()-t, gg[f], ax[f[0]], ax[f[1]], ax[f[2]]), flush=True)
    # 상위 후보 8개
    flat = gg.ravel(); order = np.argsort(flat)[:8]
    cands = []
    for o in order:
        a,b,c = np.unravel_index(o, gg.shape)
        cands.append((flat[o], ax[a], ax[b], ax[c]))
    if N == 61:
        top = cands
for gv,a,b,c in top:
    print('     후보 gap %.4e  (%.4f, %.4f, %.4f)' % (gv,a,b,c))

# ---- 3. 후보 주변 단계적 미세화 (격자만) -------------------------------------
print('\n[3] 후보 **전부** 미세화 (격자 5단계, gradient 미사용)')
for gv, a, b, c in top:          # 전부 미세화할 것. 상위 몇 개만 보면 진짜 노드를 놓친다.
    ctr = np.array([a, b, c]); half = 0.5/60
    for step in range(5):
        ax2 = [np.linspace(ctr[i]-half, ctr[i]+half, 13) for i in range(3)]
        Q = np.array(np.meshgrid(*ax2, indexing='ij')).reshape(3, -1).T
        g = gaps([list(q) for q in Q])
        j = int(np.argmin(g)); ctr = Q[j]; best = g[j]; half /= 4
    fl = freqs(ctr)
    tag = 'polar=0.5 평면으로 미끄러짐' if abs(ctr[P]-0.5) < 1e-3 else \
          ('거울면으로 미끄러짐' if any(abs(ctr[i])<1e-3 or abs(abs(ctr[i])-0.5)<1e-3 for i in mir) else '일반위치 유지')
    print('  시작 (%.4f,%.4f,%.4f) -> (%.6f,%.6f,%.6f)  gap %.4e THz  E %.5f  [%s]'
          % (a,b,c,ctr[0],ctr[1],ctr[2],best,fl[BAND-1],tag), flush=True)

# ---- 4. 기록된 점 재현 --------------------------------------------------------
print('\n[4] 기록된 점 재현 (results/weyl_trend/SUMMARY.md)')
for q in [(0.1818, 0.0563, 0.1410)]:
    g = float(gaps([list(q)])[0]); fl = freqs(q)
    print('  (%.4f,%.4f,%.4f)  gap %.4e THz  E %.5f THz  (기록: 6.67e-3)' % (*q, g, fl[BAND-1]))
