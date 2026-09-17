"""구조별 (a,b,c) 역할 -> native 축 인덱스 매핑.

c = 편극축(2_1 나사축), a = 긴 거울축, b = 짧은 거울축.
parent_pristine 의 native 순서와 일치하므로 그 구조의 seekpath 경로를 그대로 쓸 수 있고,
다른 구조에는 이 매핑으로 같은 경로를 옮긴다. POSCAR 를 건드리지 않으므로
FORCE_SETS/BORN 이 그대로 유효하다.
"""
import numpy as np


def roles(src):
    """(role2native, polar_idx, mirror_idx) 반환. role2native[i] = 역할 a,b,c 의 native 인덱스."""
    import spglib
    from phonopy.interface.calculator import read_crystal_structure
    u, _ = read_crystal_structure(filename=src + '/POSCAR', interface_mode='vasp')
    ds = spglib.get_symmetry_dataset((u.cell, u.scaled_positions, u.numbers), symprec=1e-4)
    mir = set()
    for r in ds.rotations:
        if np.count_nonzero(r - np.diag(np.diag(r))):
            continue
        d = np.diag(r)
        if np.sum(d == -1) == 1:
            mir.add(int(np.where(d == -1)[0][0]))
    p = [i for i in range(3) if i not in mir][0]
    L = np.linalg.norm(np.array(u.cell), axis=1)
    m = sorted(mir, key=lambda i: -L[i])       # 긴 거울축 먼저
    return [m[0], m[1], p], p, sorted(mir)


def to_native(k_std, r2n):
    """표준 (a,b,c) 성분을 native 성분으로."""
    out = [0.0, 0.0, 0.0]
    for role, comp in enumerate(k_std):
        out[r2n[role]] = comp
    return out


# 표준 고대칭점 (a,b,c 역할 기준)
HS = {'G': (0, 0, 0), 'X': (0.5, 0, 0), 'S': (0.5, 0.5, 0), 'Y': (0, 0.5, 0),
      'Z': (0, 0, 0.5), 'U': (0.5, 0, 0.5), 'R': (0.5, 0.5, 0.5), 'T': (0, 0.5, 0.5)}
PATH = ['G', 'X', 'S', 'Y', 'G', 'Z', 'U', 'R', 'T', 'Z']
