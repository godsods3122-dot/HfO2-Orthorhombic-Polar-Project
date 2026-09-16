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


ka, kb, Oa, Ob = field('source/parent_pristine')
mag = np.hypot(Oa, Ob)

STK = [pe.withStroke(linewidth=3.0, foreground='black')]
fig, ax = plt.subplots(figsize=(10.4, 7.2))
im = ax.pcolormesh(ka, kb, np.log10(mag + 1e-2).T, cmap='inferno',
                   shading='gouraud', rasterized=True)
ax.streamplot(ka, kb, Oa.T, Ob.T, color='white', density=1.7,
              linewidth=0.85, arrowsize=1.05)
for a, bb, c in NODES:
    ax.plot(a, bb, 'o', ms=15, mfc='none', mec='#00e5ff', mew=3.0, zorder=6)
    ax.annotate('$\\chi=%+d$' % c, xy=(a, bb),
                xytext=(0, 30 if bb > 0 else -30), textcoords='offset points',
                color='#00e5ff', ha='center', va='center', fontsize=15,
                fontweight='bold', zorder=8, path_effects=STK)
ax.set_xlim(ka[0], ka[-1]); ax.set_ylim(kb[0], kb[-1])
ax.set_xlabel('$k_a$  (reduced)')
ax.set_ylabel('$k_b$  (reduced)')
cb = fig.colorbar(im, ax=ax, pad=0.02, fraction=0.046)
cb.set_label('$\\log_{10}|\\Omega|$   ($\\AA^2$)', fontsize=13)
ax.set_title('Berry curvature field in the polar plane  ($k_c$ = 0)\n'
             'bands 1–17;  $\\Omega_c \\equiv 0$ here by $(2_1\\!\\parallel\\! c)\\cdot T$,'
             '  so the field lies in the plane', fontsize=15, pad=12)
fig.savefig('figs/fig9_berry_field.png')
print('fig9 저장.  격자 %d x %d,  |Ω| %.2e .. %.2e' % (len(ka), len(kb), mag.min(), mag.max()))
