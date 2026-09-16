#!/usr/bin/env python3
"""Fig 9: 편극 평면(k_c = 0)의 Berry curvature 벡터장.

점유 다양체(밴드 1–17)의 Berry curvature 를 Fukui–Hatsugai–Suzuki plaquette
방식으로 직접 계산한다.  Simphony 와 무관한 독립 구현이고, 같은 코드로 잰
노드 전하가 Simphony PN.out 과 일치한다 (아래 '검증').

## 왜 이 평면에서 벡터장이 되나

`(2₁∥c)·T` 가 `k_c = 0` 평면의 모든 k 를 고정시키므로 Hamiltonian 을 실수로
만들 수 있고, 그 결과 평면에 수직한 성분 `Ω_c` 가 **정확히 0** 이다.
실측: `|Ω_c| ≤ 1.3e-07` vs `|Ω_ab| ~ 8.7e+02` — 9자리 차이.
따라서 Berry curvature 가 평면 안에 누워 있고, 네 바일 노드가 그 장의
source / sink 로 나타난다.

## 검증 (이 스크립트와 같은 FHS 구현)

노드를 감싸는 직육면체 표면의 총 flux / 2π:

```
(+0.14649, +0.07085, 0)  ->  +1
(+0.14649, -0.07085, 0)  ->  -1
(-0.14649, +0.07085, 0)  ->  -1
(-0.14649, -0.07085, 0)  ->  +1
```

Simphony `WeylChirality_calc` 의 PN.out 결과와 완전히 일치한다.

⚠️ **상자를 c 방향으로 길쭉하게 잡아야 한다.**  이 노드는 c 방향으로 gap 이
2차로 매우 느리게 열려서(면내 속도의 1/350) 정육면체를 쓰면 ±c 면이 거의
축퇴 영역에 걸친다.  거기서 Berry curvature 가 폭발해 plaquette 하나가 2π 를
넘고 알리아싱이 난다 — 격자를 n=96 까지 키워도 부호가 틀린 채로 수렴한다.
면내 L 에 대해 `L_c = 5L` 로 두면 모든 크기에서 +1 로 안정된다.
(Simphony Wilson loop 이 `Nk2 = 801` 을 요구한 것과 같은 병이다.)
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from style import setup, GREY
setup()
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from weyl_scan import get_ph

NOCC = 17
NODES = [(+0.1464927, +0.0708493, +1), (+0.1464927, -0.0708493, -1),
         (-0.1464927, +0.0708493, -1), (-0.1464927, -0.0708493, +1)]


class Berry:
    """점유 다양체의 Berry curvature — FHS plaquette."""

    def __init__(self, src):
        self.ph = get_ph(src)
        A = np.array(self.ph.primitive.cell)
        self.B = 2 * np.pi * np.linalg.inv(A).T      # k_cart = k_red @ B
        self.n = np.linalg.norm(self.B, axis=1)      # |b_i|
        assert np.linalg.det(self.B / self.n[:, None]) > 0, '좌표계가 오른손이 아니다'

    def red(self, k0, d):
        return np.asarray(k0) + np.asarray(d) / self.n

    def vecs(self, qs):
        self.ph.run_qpoints(list(map(list, qs)), with_eigenvectors=True)
        return np.array(self.ph.qpoints.eigenvectors)[:, :, :NOCC]

    @staticmethod
    def link(Va, Vb):
        d = np.linalg.det(np.einsum('qij,qik->qjk', Va.conj(), Vb))
        return d / np.abs(d)

    def flux(self, corners):
        V = [self.vecs(c) for c in corners]
        u = (self.link(V[0], V[1]) * self.link(V[1], V[2])
             * self.link(V[2], V[3]) * self.link(V[3], V[0]))
        return -np.angle(u)

    def omega(self, k0, pts, h, comps=(0, 1, 2)):
        pts = np.asarray(pts, float)
        out = np.zeros_like(pts)
        for i in comps:
            j, k = (i + 1) % 3, (i + 2) % 3          # 오른손 (j,k) -> i
            ej, ek = np.eye(3)[j] * h / 2, np.eye(3)[k] * h / 2
            cs = [pts - ej - ek, pts + ej - ek, pts + ej + ek, pts - ej + ek]
            out[:, i] = self.flux([self.red(k0, c) for c in cs]) / h ** 2
        return out

    def box_chern(self, k0, L, n):
        """반변 L(축별 3성분) 직육면체 표면의 총 flux / 2π."""
        L = np.broadcast_to(np.asarray(L, float), (3,))
        tot = 0.0
        for ax in range(3):
            j, k = (ax + 1) % 3, (ax + 2) % 3
            for s in (+1, -1):
                a, b = (j, k) if s > 0 else (k, j)   # 바깥 법선 = s*e_ax
                ta = np.linspace(-L[a], L[a], n + 1)
                tb = np.linspace(-L[b], L[b], n + 1)
                P = np.zeros((n + 1, n + 1, 3))
                P[..., ax] = s * L[ax]
                P[..., a] = ta[:, None]
                P[..., b] = tb[None, :]
                cs = [P[:-1, :-1], P[1:, :-1], P[1:, 1:], P[:-1, 1:]]
                tot += self.flux([self.red(k0, c.reshape(-1, 3)) for c in cs]).sum()
        return tot / (2 * np.pi)


RINGS = 'figs/berry_rings_pp.npz'
RADII = (0.0006, 0.0025)                      # reduced
NANG = 8

CACHE = 'figs/berry_plane_pp.npz'


def field(src, na=131, nb=81, h=3e-4):
    if os.path.exists(CACHE):
        z = np.load(CACHE)
        return z['ka'], z['kb'], z['Oa'], z['Ob']
    b = Berry(src)
    ka = np.linspace(-0.26, 0.26, na)
    kb = np.linspace(-0.16, 0.16, nb)
    P = np.array([[x * b.n[0], y * b.n[1], 0.0] for x in ka for y in kb])
    O = b.omega([0.0, 0.0, 0.0], P, h, comps=(0, 1))   # Ω_c 는 대칭으로 0
    Oa, Ob = O[:, 0].reshape(na, nb), O[:, 1].reshape(na, nb)
    np.savez(CACHE, ka=ka, kb=kb, Oa=Oa, Ob=Ob)
    return ka, kb, Oa, Ob


LOCAL = 'figs/berry_local_pp.npz'


def local_field(src, nodes, r=0.0006, n=11, h=4e-5):
    """인셋용 — 노드 주변 조밀 격자 (reduced 반경 r)."""
    if os.path.exists(LOCAL):
        z = np.load(LOCAL)
        return z['ka'], z['kb'], z['Oa'], z['Ob']
    b = Berry(src)
    KA, KB, OA, OB = [], [], [], []
    for a0, b0, _ in nodes:
        xa = np.linspace(a0 - r, a0 + r, n)
        xb = np.linspace(b0 - r, b0 + r, n)
        P = np.array([[x * b.n[0], y * b.n[1], 0.0] for x in xa for y in xb])
        # ⚠️ 인셋 패치는 monopole 영역 *안* 이라 패치 평균이 배경이 아니다.
        # 빼면 monopole 자체를 깎아 균일장처럼 만들어 버린다.  여기서는 안 뺀다.
        # (링은 원 둘레 평균이 정확히 배경이므로 거기서만 뺀다.)
        O = b.omega([0.0, 0.0, 0.0], P, h, comps=(0, 1))
        KA.append(xa); KB.append(xb)
        OA.append(O[:, 0].reshape(n, n)); OB.append(O[:, 1].reshape(n, n))
    np.savez(LOCAL, ka=np.array(KA), kb=np.array(KB),
             Oa=np.array(OA), Ob=np.array(OB))
    return np.array(KA), np.array(KB), np.array(OA), np.array(OB)


# 부호 규약.  이 FHS 구현의 Ω 부호는 Simphony 의 χ 와 전체적으로 반대다
# (프로젝트 노트의 "up to a global sign" 그대로).  근거리 측정:
#   (+a,+b) <cos>=-1.0000  vs Simphony χ=+1
#   (+a,-b) <cos>=+1.0000  vs Simphony χ=-1
#   (-a,+b) <cos>=+1.0000  vs Simphony χ=-1
#   (-a,-b) <cos>=-1.0000  vs Simphony χ=+1
# 교대 패턴은 완벽히 일치하므로 어느 노드가 source/sink 인지는 이 계산이
# 독립적으로 확인한다.  전체 부호만 Simphony 에 맞춰 뒤집는다.
SIGN = -1.0

def ring_field(src, nodes):
    """노드마다 로그 간격 링 위의 Ω.

    주 격자 간격은 0.004 인데 monopole 영역은 반경 0.0008 이라 격자점이 그
    안에 하나도 안 들어간다.  그래서 큰 그림에서는 노드가 안장점처럼 보인다.
    링을 따로 계산해 얹으면 노드 바로 옆의 방사형 구조가 드러난다.
    """
    if os.path.exists(RINGS):
        z = np.load(RINGS)
        return z['P'], z['O']
    b = Berry(src)
    P, O = [], []
    th = np.linspace(0, 2 * np.pi, NANG, endpoint=False)
    for a0, b0, _ in nodes:
        for r in RADII:
            d = np.c_[r * np.cos(th) * b.n[0], r * np.sin(th) * b.n[1],
                      np.zeros(NANG)]
            o = b.omega([a0, b0, 0.0], d, min(4e-5, r * b.n[0] * 0.12),
                        comps=(0, 1))
            # 국소 배경 제거.  monopole 성분은 원 둘레 평균이 정확히 0 이므로
            # 링 평균이 곧 배경이다.  빼고 나면 r=0.0006 에서 화살표가 100 %
            # 바깥(또는 안)을 향한다 (빼기 전 79 %).
            o[:, :2] -= o[:, :2].mean(0)
            P.append(np.c_[a0 + r * np.cos(th), b0 + r * np.sin(th)])
            O.append(o[:, :2])
    P, O = np.vstack(P), np.vstack(O)
    np.savez(RINGS, P=P, O=O)
    return P, O


ka, kb, Oa, Ob = field('source/parent_pristine')
rP, rO = ring_field('source/parent_pristine', NODES)
rO = SIGN * rO
Oa, Ob = SIGN * Oa, SIGN * Ob
INS = [NODES[0], NODES[1]]                      # χ=+1 하나, χ=−1 하나
lka, lkb, lOa, lOb = local_field('source/parent_pristine', INS)
lOa, lOb = SIGN * lOa, SIGN * lOb

RED, BLU = '#d62728', '#1f77b4'
ARROW = '#3b8fd4'


def uquiver(ax, A, B, u, v, lo=8, hi=99, floor=0.18, **kw):
    """길이를 log|Ω| 로 압축해 넣은 화살표.

    |Ω| 는 노드에서 1/r² 로 발산하고 먼 곳에서는 거의 0 이라 선형 길이로는
    쓸 수 없다.  log 를 백분위로 정규화해 [floor, 1] 로 눌러 넣는다.
    꼬리를 길게 보이도록 머리를 작게 잡았다.
    """
    m = np.hypot(u, v)
    good = m > 0
    L = np.full_like(m, floor)
    if good.any():
        lg = np.log10(m[good])
        a, b = np.percentile(lg, [lo, hi])
        L[good] = np.clip((lg - a) / max(b - a, 1e-12), 0.0, 1.0) * (1 - floor) + floor
    mm = np.where(good, m, 1.0)
    ax.quiver(A, B, u / mm * L, v / mm * L, color=ARROW, angles='xy',
              scale_units='width', headwidth=4.0, headlength=4.0,
              headaxislength=3.4, **kw)


fig, ax = plt.subplots(figsize=(11.0, 7.6))
st = 3
A, B = np.meshgrid(ka[::st], kb[::st], indexing='ij')
U, V = Oa[::st, ::st].copy(), Ob[::st, ::st].copy()
near = np.zeros(A.shape, bool)                 # 링이 덮는 자리는 격자를 뺀다
for a0, b0, _ in NODES:
    near |= np.hypot(A - a0, B - b0) < 0.006
U[near] = 0.0; V[near] = 0.0
uquiver(ax, A, B, U, V, scale=30, width=0.0022, zorder=2)
uquiver(ax, rP[:, 0], rP[:, 1], rO[:, 0].copy(), rO[:, 1].copy(),
        scale=30, width=0.0026, zorder=3, lo=2, hi=98, floor=0.35)

LBL = {0: (34, 14), 1: (34, -14), 2: (0, 30), 3: (0, -30)}
for n, (a, bb, c) in enumerate(NODES):
    ax.plot(a, bb, 'o', ms=14, color=RED if c > 0 else BLU,
            mec='white', mew=1.5, zorder=7)
    ax.annotate('$W_%d$' % (n + 1), xy=(a, bb), xytext=LBL[n],
                textcoords='offset points', color='#111111', ha='center',
                va='center', fontsize=17, fontweight='bold', zorder=8)

# 확대 인셋 — 노드가 없는 가운데 세로 띠에 두고 연결선으로 잇는다
for (a, bb, c), xa, xb, ua, ub, loc in zip(
        INS, lka, lkb, lOa, lOb,
        ([0.385, 0.605, 0.225, 0.285], [0.385, 0.090, 0.225, 0.285])):
    axi = ax.inset_axes(loc, zorder=12)
    axi.set_facecolor('white')
    axi.patch.set_alpha(1.0)
    A2, B2 = np.meshgrid(xa, xb, indexing='ij')
    uquiver(axi, A2, B2, ua.copy(), ub.copy(), scale=13, width=0.0085,
            lo=5, hi=95, floor=0.48)
    axi.plot(a, bb, 'o', ms=15, color=RED if c > 0 else BLU,
             mec='white', mew=2.0, zorder=6)
    axi.set_xlim(xa[0], xa[-1]); axi.set_ylim(xb[0], xb[-1])
    axi.set_xticks([]); axi.set_yticks([])
    for sp in axi.spines.values():
        sp.set_color('#555555'); sp.set_linewidth(1.4)
    ax.indicate_inset_zoom(axi, edgecolor='#555555', alpha=0.95, lw=1.3)

ax.set_xlim(ka[0], ka[-1]); ax.set_ylim(kb[0], kb[-1])
ax.set_xlabel('$k_a$  (reduced)')
ax.set_ylabel('$k_b$  (reduced)')
fig.savefig('figs/fig9_berry_field.png')
print('fig9 저장.  주 격자 %d x %d,  인셋 %d x %d' % (len(ka), len(kb), lka.shape[1], lkb.shape[1]))
