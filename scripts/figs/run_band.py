#!/usr/bin/env python3
"""Simphony BulkBand 을 표준 경로(c=편극축)로 돌리고 bulkek.dat 을 남긴다."""
import argparse, os, shutil, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from axes import roles, to_native, HS, PATH

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def build(workdir, src, nk, segs=None):
    r2n, _, _ = roles(src)
    segs = segs or [(PATH[i], PATH[i + 1]) for i in range(len(PATH) - 1)]
    lines = ['KPATH_BULK', str(len(segs))]
    for A, B in segs:
        a = to_native(HS[A], r2n); b = to_native(HS[B], r2n)
        lines.append('%s %14.8f %14.8f %14.8f   %s %14.8f %14.8f %14.8f'
                     % (A, a[0], a[1], a[2], B, b[0], b[1], b[2]))
    p = os.path.join(workdir, 'pn.in')
    L = open(p).read().splitlines(); out = []; i = 0
    while i < len(L):
        x = L[i]; t = x.strip()
        if t == 'KPATH_BULK':
            out += lines
            i += 2 + int(L[i + 1]); continue
        if t.startswith(('BulkBand_calc', 'BulkGap_cube_calc', 'WeylChirality_calc',
                         'SlabBand_calc', 'SlabArc_calc', 'SlabSS_calc')):
            key = t.split('=')[0].strip()
            out.append('  %s = %s' % (key, 'T' if key == 'BulkBand_calc' else 'F')); i += 1; continue
        if t.startswith('NumOccupied'): out.append('  NumOccupied = 17'); i += 1; continue
        if 'Nk1 = ' in x: out.append('  Nk1 = %d' % nk); i += 1; continue
        if 'LOTO_method' in x: i += 1; continue
        if 'Gap_threshold' in x:
            out.append("  LOTO_method = 'phonopy'")
            out.append('  Gap_threshold = 0.00001'); i += 1; continue
        out.append(x); i += 1
    open(p, 'w').write('\n'.join(out) + '\n')
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True); ap.add_argument('--src', required=True)
    ap.add_argument('--nk', type=int, default=120)
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    segs = build(a.dir, a.src, a.nk)
    print('경로 %d 구간, 구간당 %d 점' % (len(segs), a.nk))
    for f in ('bulkek.dat',):
        q = os.path.join(a.dir, f)
        if os.path.exists(q): os.remove(q)
    subprocess.run([os.path.join(ROOT, 'simphony', 'src', 'pn.x')], cwd=a.dir,
                   stdout=open(os.path.join(a.dir, 'run.log'), 'w'),
                   stderr=subprocess.STDOUT, check=False)
    src_dat = os.path.join(a.dir, 'bulkek.dat')
    if not os.path.exists(src_dat):
        print('bulkek.dat 없음 - run.log 확인'); return
    if a.out:
        shutil.copy(src_dat, a.out); print('저장:', a.out)
    d = np.loadtxt(src_dat)
    print('점 %d, 밴드 %d' % (len(d) // 36, 36))


if __name__ == '__main__':
    main()
