#!/usr/bin/env python3
"""Weyl phonon chirality — 구면 Berry flux (Fukui-Hatsugai-Suzuki).

Simphony 없이 phonopy 고유벡터만으로 χ 를 낸다. parent_pristine 에서 Simphony 와
개별값·상대부호·합이 일치함을 확인한 구현이다 (전체 부호 규약만 반대).

판정은 gap 크기가 아니라 **χ 가 여러 반지름에서 정수로 안정한가**로 한다.
대칭이 축퇴를 강제하는 자리(nodal plane, 거울선)에서는 구가 축퇴면에 관통되어
χ 자체가 정의되지 않으므로, 거기서 나온 숫자는 인용하면 안 된다.

사용: python3 weyl_chirality_flux.py SRC BAND K1 K2 K3 [R ...]
      BAND=17 이면 1..17 을 occupied 로 보고 17-18 사이 gap 의 χ 를 낸다.
"""
import sys
import numpy as np

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weyl_scan import get_ph

SRC = sys.argv[1]
NOCC = int(sys.argv[2])
CEN = np.array([float(x) for x in sys.argv[3:6]])
RADII = [float(x) for x in sys.argv[6:]] or [0.002, 0.004, 0.008]

ph = get_ph(SRC)
rec = np.linalg.inv(ph.primitive.cell).T * 2 * np.pi     # rows = b_i, 1/Angstrom


def flux(center, r_cart, nth=40, nph=80):
    th = np.linspace(0, np.pi, nth)
    phi = np.linspace(0, 2 * np.pi, nph, endpoint=False)
    T, P = np.meshgrid(th, phi, indexing='ij')
    dirs = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
    kred = (dirs * r_cart) @ np.linalg.inv(rec)
    Q = (center + kred).reshape(-1, 3)
    ph.run_qpoints([list(q) for q in Q], with_eigenvectors=True)
    V = ph.qpoints.eigenvectors
    nb = V.shape[1]
    V = V.reshape(nth, nph, nb, nb)[:, :, :, :NOCC]      # occupied block

    def link(A, B):
        M = np.einsum('...ji,...jk->...ik', A.conj(), B)
        d = np.linalg.det(M)
        return d / np.abs(d)

    Vp = np.concatenate([V, V[:, :1]], axis=1)           # phi 주기 닫기
    U1 = link(Vp[:-1, :-1], Vp[:-1, 1:])
    U2 = link(Vp[:-1, :-1], Vp[1:, :-1])
    U3 = link(Vp[1:, :-1], Vp[1:, 1:])
    U4 = link(Vp[:-1, 1:], Vp[1:, 1:])
    return np.angle(U1 * U4 / (U3 * U2)).sum() / (2 * np.pi)


ph.run_qpoints([list(CEN)])
f = ph.qpoints.frequencies[0]
print('k = (%.7f, %.7f, %.7f)' % tuple(CEN))
print('gap(%d-%d) = %.4e THz    E = %.6f THz' % (NOCC, NOCC + 1, f[NOCC] - f[NOCC - 1], f[NOCC - 1]))
for r in RADII:
    print('  r = %.4f 1/Ang   chi = %+.4f' % (r, flux(CEN, r)))
