"""검증 (C): 구면 Wilson loop / WCC 감김수 — FHS plaquette flux 와 다른 알고리즘.

구를 극각 theta 로 자르고, 각 위도원에서 점유 17밴드의 Wilson loop 를 만들어
고유값 위상(=WCC)을 뽑는다. theta: 0 -> pi 로 갈 때 WCC 합의 총 감김수가
Chern 수 = chirality 다. Berry flux 합산과 수학적으로 동등하지만 구현 경로가 다르고,
게이지 선택·플라케트 부호 규약에 의존하지 않는다.
"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/',2)[0])
from weyl_scan import get_ph

SRC=sys.argv[1]; K0=np.array([float(x) for x in sys.argv[2:5]])
NOCC=17
R=float(sys.argv[5]) if len(sys.argv)>5 else 0.004
NTH=int(sys.argv[6]) if len(sys.argv)>6 else 61
NPH=int(sys.argv[7]) if len(sys.argv)>7 else 80
ph=get_ph(SRC)
rec=np.linalg.inv(ph.primitive.cell).T*2*np.pi

def occ(Q):
    ph.run_qpoints([list(q) for q in Q],with_eigenvectors=True)
    return ph.qpoints.eigenvectors[:,:,:NOCC]

tot=0.0; prev=None; wind=0.0
phis=np.linspace(0,2*np.pi,NPH,endpoint=False)
ths=np.linspace(1e-3,np.pi-1e-3,NTH)
sums=[]
for th in ths:
    d=np.stack([np.sin(th)*np.cos(phis),np.sin(th)*np.sin(phis),np.cos(th)*np.ones(NPH)],-1)
    Q=K0+(d*R)@np.linalg.inv(rec)
    V=occ(Q)
    W=np.eye(NOCC,dtype=complex)
    for i in range(NPH):
        j=(i+1)%NPH
        W=W@(V[i].conj().T@V[j])
    ev=np.linalg.eigvals(W)
    s=np.angle(np.linalg.det(W))          # WCC 합 = det 의 위상
    sums.append(s)
sums=np.array(sums)
# 감김수: 인접 theta 사이 위상차를 (-pi,pi] 로 펴서 누적
d=np.diff(sums); d=(d+np.pi)%(2*np.pi)-np.pi
chern=d.sum()/(2*np.pi)
print('k0=(%.7f,%.7f,%.7f)  R=%.4f 1/Ang  (NTH=%d, NPH=%d)'%(*K0,R,NTH,NPH))
print('  WCC 합의 총 감김수 = %+.4f   -> chi = %+d'%(chern,int(round(chern))))
