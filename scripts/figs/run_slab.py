#!/usr/bin/env python3
"""Simphony 표면 계산 실행: SlabSS_calc (경로별 표면 스펙트럼) 또는 SlabArc_calc (등에너지 아크맵).

SURFACE 는 (1 0 0)/(0 1 0) 로 두어 법선이 c=편극축이 되게 한다 — 이 Weyl 계열의
아크가 나오는 유일한 배향이다 (results/domain_route/SUMMARY.md D절).
"""
import argparse, os, subprocess, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def patch(path, mode, nk1, nk2, omin, omax, onum, earc, eta, npl, kpath=None, kplane=None):
    L = open(path).read().splitlines(); out = []; i = 0
    while i < len(L):
        x = L[i]; t = x.strip()
        if t == 'KPLANE_SLAB' and kplane:
            out += ['KPLANE_SLAB'] + list(kplane)
            i += 4; continue
        if t == 'KPATH_SLAB' and kpath:
            out += ['KPATH_SLAB', str(len(kpath))] + kpath
            i += 2 + int(L[i + 1]); continue
        if t.startswith(('BulkBand_calc', 'BulkGap_cube_calc', 'WeylChirality_calc',
                         'SlabBand_calc', 'SlabArc_calc', 'SlabSS_calc')):
            key = t.split('=')[0].strip()
            on = (key == 'SlabSS_calc' and mode == 'ss') or (key == 'SlabArc_calc' and mode == 'arc')
            out.append('  %s = %s' % (key, 'T' if on else 'F')); i += 1; continue
        if t.startswith('NumOccupied'): out.append('  NumOccupied = 17'); i += 1; continue
        if t.startswith('NSLAB'): out.append('  NSLAB = 20'); i += 1; continue
        if 'Nk1 = ' in x: out.append('  Nk1 = %d' % nk1); i += 1; continue
        if 'Nk2 = ' in x: out.append('  Nk2 = %d' % nk2); i += 1; continue
        if 'NP = ' in x: out.append('  NP = %d' % npl); i += 1; continue
        if 'Eta_Arc' in x: out.append('  Eta_Arc = %g' % eta); i += 1; continue
        if 'E_arc' in x: out.append('  E_arc = %.6f' % earc); i += 1; continue
        if 'OmegaNum' in x: out.append('  OmegaNum = %d' % onum); i += 1; continue
        if 'OmegaMin' in x: out.append('  OmegaMin = %.5f' % omin); i += 1; continue
        if 'OmegaMax' in x: out.append('  OmegaMax = %.5f' % omax); i += 1; continue
        if 'LOTO_method' in x: i += 1; continue
        if 'Gap_threshold' in x:
            out.append("  LOTO_method = 'phonopy'")
            out.append('  Gap_threshold = 0.00001'); i += 1; continue
        out.append(x); i += 1
    open(path, 'w').write('\n'.join(out) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dir', required=True)
    ap.add_argument('--mode', choices=['ss', 'arc'], required=True)
    ap.add_argument('--nk1', type=int, default=160); ap.add_argument('--nk2', type=int, default=160)
    ap.add_argument('--omin', type=float, default=9.60); ap.add_argument('--omax', type=float, default=10.60)
    ap.add_argument('--onum', type=int, default=300)
    ap.add_argument('--earc', type=float, default=10.0869)
    ap.add_argument('--eta', type=float, default=0.004)
    ap.add_argument('--np', dest='npl', type=int, default=5)
    ap.add_argument('--kpath', nargs='*', default=None)
    ap.add_argument('--kplane', nargs='*', default=None)
    a = ap.parse_args()
    patch(os.path.join(a.dir, 'pn.in'), a.mode, a.nk1, a.nk2, a.omin, a.omax,
          a.onum, a.earc, a.eta, a.npl, a.kpath, a.kplane)
    subprocess.run([os.path.join(ROOT, 'simphony', 'src', 'pn.x')], cwd=a.dir,
                   stdout=open(os.path.join(a.dir, 'run.log'), 'w'),
                   stderr=subprocess.STDOUT, check=False)
    for f in sorted(os.listdir(a.dir)):
        if f.startswith(('dos.dat', 'arc.dat')):
            print('  ', f, os.path.getsize(os.path.join(a.dir, f)))


if __name__ == '__main__':
    main()
