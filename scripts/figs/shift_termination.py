#!/usr/bin/env python3
"""표면 종단면을 바꾸기 위해 원자를 편극축 방향으로 통째로 민다.

WannierTools/Simphony 는 슬랩을 단위포 경계에서 자른다.  원자를 면직 방향으로
δ 만큼 **강체 이동**시키면 어느 면이 잘리는지가 바뀐다.  H(R) 은 궤도
i(τ_i+δ+R) 와 j(τ_j+δ) 사이의 값이라 강체 이동에 불변이고, LO-TO 쌍극자항도
τ 차이에만 의존하므로 hr.dat / BORN 을 다시 만들 필요가 없다.

⚠️ Simphony 가 실제로 읽는 좌표는 **`pn.in` 의 `ATOMS` 카드**다.  POSCAR 만
고치면 아무 일도 일어나지 않는다 (이걸로 한 번 헛돌았다).

사용:
    python3 shift_termination.py <pn.in> <axis 0|1|2> <shift>      # 제자리 수정
"""
import sys


def main():
    path, axis, sh = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    L = open(path).read().splitlines()
    try:
        i = next(k for k, x in enumerate(L) if x.strip() == 'ATOMS')
    except StopIteration:
        sys.exit('ATOMS 카드를 찾을 수 없다: %s' % path)
    nat = int(L[i + 1].split()[0])
    mode = L[i + 2].strip().lower()
    if not mode.startswith(('d', 'direct')):
        sys.exit('Direct 좌표계만 지원한다 (현재: %r)' % L[i + 2])
    zs = []
    for k in range(i + 3, i + 3 + nat):
        p = L[k].split()
        v = [float(p[2]), float(p[3]), float(p[4])]
        v[axis] = (v[axis] + sh) % 1.0
        zs.append(v[axis])
        L[k] = '%-8s %12.5f %24.10f %24.10f %24.10f' % (p[0], float(p[1]), *v)
    open(path, 'w').write('\n'.join(L) + '\n')
    u = sorted(set(round(z, 4) for z in zs))
    print('%s   axis %d, shift %+.4f' % (path, axis, sh))
    print('  이동 후 층:', u)
    g = [round(u[(j + 1) % len(u)] + (1 if j == len(u) - 1 else 0) - u[j], 4)
         for j in range(len(u))]
    print('  빈틈:', g, '  경계(z=0)를 가로지르는 빈틈 =', g[-1])


if __name__ == '__main__':
    main()
