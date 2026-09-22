#!/usr/bin/env python3
"""m1_mirror (미러 −1 %) 기준 그림 데이터 일괄 생성 — phonopy 직접 계산.

지금까지 fig2/fig5/fig9 는 parent_pristine(일반 구조), fig8 은 미러 계열을 써서
두 계열이 섞여 있었다. 두 데이터셋은 주파수가 0.14 THz 어긋나 비교 대상이 아니다
(results/strain_audit/PARENT_VS_MIRROR.md). 전부 미러 −1 % 로 통일한다.

출력: figs/m1_figdata.npz
  path_x, path_E, ticks, labels     표준 경로 Γ-X-S-Y-Γ-Z-U-R-T-Z 밴드
  gamma_break                       안쪽 Γ 인덱스 (LO-TO 방향 의존 불연속)
  xc_t, xc_E                        Γ-X 위 band17/18 교점
  cut_k_a, cut_E_a / cut_k_b, cut_E_b   노드 통과 1D 절단 (절대 좌표)
  cone_a, cone_E17, cone_E18        노드 주변 (k_a,k_b) 격자
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from axes import roles, to_native, HS, PATH
from weyl_scan import get_ph, gap_at

SRC = 'source/m1_mirror'
BAND = 17
NK = 120                       # 구간당 점수
NODE_STD = (0.1294182, 0.0637717, 0.0)     # 표준 (a,b,c) — plane_survey.py 로 재확인
HALF = 0.035
CONE_R, CONE_N = 0.0085, 101

r2n, pol, mir = roles(SRC)
ph = get_ph(SRC)
print('m1_mirror  role2native=%s  polar native=%d' % (r2n, pol))

def native(kstd):
    return to_native(kstd, r2n)

def freqs(qs, chunk=20000):
    out = []
    for i in range(0, len(qs), chunk):
        ph.run_qpoints([list(q) for q in qs[i:i+chunk]])
        out.append(ph.qpoints.frequencies)
    return np.concatenate(out)

# ---- 1) 표준 경로 --------------------------------------------------------------
rec = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi       # native 기준 역격자
segs, xs, x0 = [], [], 0.0
Q = []
for a, b in zip(PATH[:-1], PATH[1:]):
    ka, kb = np.array(HS[a], float), np.array(HS[b], float)
    L = np.linalg.norm((np.array(native(kb)) - np.array(native(ka))) @ rec)
    t = np.linspace(0, 1, NK, endpoint=False)
    for tt in t:
        Q.append(native(ka + tt * (kb - ka)))
    xs.append(x0 + t * L); segs.append(x0); x0 += L
Q.append(native(np.array(HS[PATH[-1]], float)))
xs = np.concatenate(xs); xs = np.append(xs, x0)
F = freqs(Q)
ticks = np.array(segs + [x0])
print('경로 %d점, 길이 %.4f' % (len(Q), x0))

# ---- 2) Γ-X 위 교점 ------------------------------------------------------------
ts = np.linspace(0, 0.5, 2001)
g = gap_at(ph, [list(native((t, 0, 0))) for t in ts], BAND)
i = int(np.argmin(g[:1600])); lo, hi = ts[i] - 0.001, ts[i] + 0.001
for _ in range(6):
    tt = np.linspace(lo, hi, 41)
    gg = gap_at(ph, [list(native((t, 0, 0))) for t in tt], BAND)
    j = int(np.argmin(gg)); lo, hi = tt[max(j-1, 0)], tt[min(j+1, 40)]
xc_t = tt[j]
xc_E = freqs([native((xc_t, 0, 0))])[0][BAND-1]
print('Γ-X 교점  t=%.6f  gap=%.3e  E=%.6f THz' % (xc_t, gg[j], xc_E))

# ---- 3) 노드 통과 1D 절단 (절대 좌표) -------------------------------------------
cuts = {}
for ax_i, name in ((0, 'a'), (1, 'b')):
    k = np.linspace(NODE_STD[ax_i] - HALF, NODE_STD[ax_i] + HALF, 241)
    QQ = []
    for v in k:
        s = list(NODE_STD); s[ax_i] = v
        QQ.append(native(s))
    cuts['cut_k_' + name] = k
    cuts['cut_E_' + name] = freqs(QQ)
    gmin = (cuts['cut_E_' + name][:, BAND] - cuts['cut_E_' + name][:, BAND-1]).min()
    print('절단 k_%s  최소 gap %.3e' % (name, gmin))

# ---- 4) 콘 -------------------------------------------------------------------
a = np.linspace(-CONE_R, CONE_R, CONE_N)
QQ = [native((NODE_STD[0] + dx, NODE_STD[1] + dy, 0.0)) for dx in a for dy in a]
Fc = freqs(QQ)
E17 = Fc[:, BAND-1].reshape(CONE_N, CONE_N)
E18 = Fc[:, BAND].reshape(CONE_N, CONE_N)
print('콘 최소 gap %.3e  E=%.6f' % ((E18 - E17).min(), E17.ravel()[np.argmin(E18 - E17)]))

np.savez('figs/m1_figdata.npz', path_x=xs, path_E=F.T, ticks=ticks,
         labels=np.array(PATH), gamma_break=4 * NK, xc_t=xc_t, xc_E=xc_E,
         node=np.array(NODE_STD), cone_a=a, cone_E17=E17, cone_E18=E18, **cuts)
print('saved figs/m1_figdata.npz')
