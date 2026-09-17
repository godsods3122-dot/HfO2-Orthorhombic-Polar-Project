#!/usr/bin/env python3
"""hr.dat 에서 슬랩을 직접 만들어 ASR 위반을 정량화하고 교정한다.

Simphony 가 하는 것과 같은 방식(법선 방향 R 성분으로 자르기)으로 슬랩을 만든 뒤
- 교정 전: 표면층에서 sum_j Phi_ij != 0  -> Gamma 에서 음향모드가 0 이 아님
- 교정 후: 대각 자기항을 sum_j 로 다시 맞춰 ASR 복원
"""
import numpy as np


def read_hr(path):
    L = open(path).read().splitlines()
    nw = int(L[1]); nr = int(L[2])
    deg = []
    i = 3
    while len(deg) < nr:
        deg += [int(x) for x in L[i].split()]; i += 1
    deg = np.array(deg[:nr], float)
    R = np.zeros((nr, 3), int)
    H = np.zeros((nr, nw, nw), complex)
    for ir in range(nr):
        for a in range(nw * nw):
            p = L[i].split(); i += 1
            if a == 0:
                R[ir] = [int(p[0]), int(p[1]), int(p[2])]
            m = int(p[3]) - 1; n = int(p[4]) - 1
            H[ir, m, n] = float(p[5]) + 1j * float(p[6])
    return R, H, deg


def slab_blocks(R, H, deg, normal):
    """법선 방향 R 성분별로 묶는다: blocks[d] = sum_{R: R[normal]=d} H(R)/deg (면내 k=0)."""
    out = {}
    for ir in range(len(R)):
        d = int(R[ir, normal])
        out.setdefault(d, np.zeros_like(H[0]))
        out[d] += H[ir] / deg[ir]
    return out


def build_slab(blocks, N):
    nw = next(iter(blocks.values())).shape[0]
    D = np.zeros((N * nw, N * nw), complex)
    for n in range(N):
        for d, B in blocks.items():
            m = n + d
            if 0 <= m < N:
                D[n * nw:(n + 1) * nw, m * nw:(m + 1) * nw] += B
    return D


def asr_fix(D, N, nw, sqm):
    """질량가중 D 에 대한 ASR: sum_j sqrt(m_j) D_ij sqrt(m_j) 형태로 되돌려 보정."""
    D = D.copy()
    M = np.tile(sqm, N)                       # sqrt(mass) per orbital
    Phi = D * np.outer(M, M)                  # 힘상수로 환원
    for i in range(N * nw):
        blk = i // 3 * 3
        pass
    # 원자(3궤도) 단위로 행 합이 0 이 되도록 대각 3x3 블록 보정
    na = N * nw // 3
    for ia in range(na):
        s = slice(3 * ia, 3 * ia + 3)
        rowsum = Phi[s, :].reshape(3, na, 3).sum(axis=1)
        Phi[s, s] -= rowsum
    return Phi / np.outer(M, M)


def freqs(D):
    w = np.linalg.eigvalsh((D + D.conj().T) / 2)
    return np.sign(w) * np.sqrt(np.abs(w))


if __name__ == '__main__':
    import sys
    from phonopy.interface.calculator import read_crystal_structure
    src, hr, normal, N = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    u, _ = read_crystal_structure(filename=src + '/POSCAR', interface_mode='vasp')
    masses = np.array(u.masses)
    sqm = np.repeat(np.sqrt(masses), 3)

    R, H, deg = read_hr(hr)
    nw = H.shape[1]
    blocks = slab_blocks(R, H, deg, normal)
    print('법선 축 index=%d, R 성분 범위 = %s' % (normal, sorted(blocks)))

    D = build_slab(blocks, N)
    f = freqs(D)
    print('\n[교정 전] NSLAB=%d, 면내 k=0 최저 6개 모드 (THz):' % N)
    print('   ', np.round(np.sort(f)[:6], 5))

    Dc = asr_fix(D, N, nw, sqm)
    fc = freqs(Dc)
    print('\n[ASR 교정 후] 최저 6개 모드 (THz):')
    print('   ', np.round(np.sort(fc)[:6], 5))

    # 벌크 참조: 모든 층이 같으므로 주기적 합
    Bsum = sum(blocks.values())
    print('\n[벌크 참조] 주기적 Gamma 최저 6개 (THz):')
    print('   ', np.round(np.sort(freqs(Bsum))[:6], 5))
