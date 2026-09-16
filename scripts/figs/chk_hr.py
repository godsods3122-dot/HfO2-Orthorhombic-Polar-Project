"""검증 (A): phonopy 동역학행렬을 전혀 쓰지 않고, hr.dat 에서 numpy 로 직접 D(k).

hr_total.dat 은 phonopy2TBDAT.py 가 쓴 질량가중 힘상수 x factor^2 이므로
  D(k) = sum_R H(R)/deg * exp(2 pi i k.R),  sqrt(eig) = THz
파일 형식도 코드 경로도 phonopy 동역학행렬과 완전히 다르다.
NAC 없는 phonopy 와 비교해 FORCE_SETS -> D(k) 경로를 독립 검증한다.
"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/',2)[0])
sys.path.insert(0, __file__.rsplit('/',1)[0])
from slab_asr import read_hr
from phonopy import Phonopy
from phonopy.file_IO import parse_FORCE_SETS
from phonopy.interface.calculator import read_crystal_structure

SRC='/home/user/HfO2-Orthorhombic-Polar-Project/source/m1_mirror'
HR ='/home/user/HfO2-Orthorhombic-Polar-Project/work/m1_mirror/hr_total.dat'
R,H,deg=read_hr(HR); Hn=H/deg[:,None,None]

def omega_numpy(k):
    D=np.tensordot(np.exp(2j*np.pi*(R@np.asarray(k,float))),Hn,axes=(0,0))
    D=(D+D.conj().T)/2
    w=np.linalg.eigvalsh(D)
    return np.sign(w)*np.sqrt(np.abs(w))

u,_=read_crystal_structure(filename=SRC+'/POSCAR',interface_mode='vasp')
ph=Phonopy(u,supercell_matrix=np.diag((2,2,2)),primitive_matrix='P')
ph.dataset=parse_FORCE_SETS(filename=SRC+'/FORCE_SETS')
ph.produce_force_constants(calculate_full_force_constants=True)
ph.symmetrize_force_constants_by_space_group(); ph.symmetrize_force_constants(level=3)
# NAC 없음

rng=np.random.default_rng(3)
ks=[rng.random(3)-0.5 for _ in range(12)]+[[0.0637716,0.0,0.1294181],[0.4189473,0.3952085,0.2625841]]
ph.run_qpoints([list(k) for k in ks]); F=ph.qpoints.frequencies
d=max(np.abs(omega_numpy(k)-F[i]).max() for i,k in enumerate(ks))
print('[A] numpy(hr_total.dat)  vs  phonopy(NAC 없음),  k 14점 x 36밴드')
print('    최대 |차이| = %.3e THz   -> %s'%(d,'일치' if d<1e-6 else '불일치'))

print('\n[B] NAC 를 끄면 노드가 어떻게 되나 (노드가 NAC 인공물인지 확인)')
for lab,k0 in [('polar=0 계열',np.array([0.0637716,0.0,0.1294181])),
               ('일반위치 계열',np.array([0.4189473,0.3952085,0.2625841]))]:
    g0=lambda k: (lambda f: f[17]-f[16])(omega_numpy(k))
    print('  %s  기존 노드 자리에서 NAC 없는 gap = %.4e THz'%(lab,g0(k0)))
    ctr=k0.copy(); half=0.03
    for _ in range(7):
        ax=[np.linspace(ctr[i]-half,ctr[i]+half,11) for i in range(3)]
        Q=np.array(np.meshgrid(*ax,indexing='ij')).reshape(3,-1).T
        gv=np.array([g0(q) for q in Q]); j=int(np.argmin(gv)); ctr=Q[j]; best=gv[j]; half/=3
    print('     -> 근처 최소 (%.6f, %.6f, %.6f)  gap %.4e  이동량 %.4f rlu'
          %(*ctr,best,np.linalg.norm(ctr-k0)))
