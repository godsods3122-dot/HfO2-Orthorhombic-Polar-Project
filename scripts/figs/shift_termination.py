#!/usr/bin/env python3
"""표면 종단면을 바꾸기 위해 POSCAR 의 원자를 편극축 방향으로 통째로 민다.

WannierTools/Simphony 는 슬랩을 단위포 경계에서 자른다.  원자를 c 방향으로
δ 만큼 **강체 이동**시키면 어느 면이 잘리는지가 바뀐다.  H(R) 은 궤도 i(τ_i+δ+R)
와 j(τ_j+δ) 사이의 값이라 강체 이동에 불변이고, LO-TO 쌍극자항도 τ 차이에만
의존하므로 hr.dat / BORN 을 다시 만들 필요가 없다.

사용:
    python3 shift_termination.py <src POSCAR> <dst POSCAR> <axis 0|1|2> <shift>
"""
import sys


def main():
    src, dst, axis, sh = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4])
    L = open(src).read().splitlines()
    nat = sum(int(x) for x in L[6].split())
    if not L[7].strip().lower().startswith(('d', 'direct')):
        sys.exit('Direct 좌표계 POSCAR 만 지원한다 (현재: %r)' % L[7])
    out = L[:8]
    for i in range(8, 8 + nat):
        p = L[i].split()
        v = [float(p[0]), float(p[1]), float(p[2])]
        v[axis] = (v[axis] + sh) % 1.0
        out.append('  %.16f  %.16f  %.16f' % tuple(v))
    out += L[8 + nat:]
    open(dst, 'w').write('\n'.join(out) + '\n')
    z = sorted(round((float(L[i].split()[axis]) + sh) % 1.0, 4) for i in range(8, 8 + nat))
    print('%s -> %s   axis %d, shift %+.4f' % (src, dst, axis, sh))
    print('  이동 후 층 z:', sorted(set(z)))


if __name__ == '__main__':
    main()
