"""chirality 가 정말 안정한가 — 반지름 6개 x 방법 2개.

parent_pristine 에서 flux 가 r=0.005 와 0.008 사이에서 뒤집히고 WCC 는 0 을 주는
것을 발견해, pristine/m1/p1 세 구조에 대해 같은 검사를 한다.
안정하지 않으면 그림에 chi 를 적으면 안 된다.
"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0])
from weyl_scan import get_ph, gap_at
NOCC=17
def flux(ph,rec,c,r,nth=60,nph=120):
    th=np.linspace(0,np.pi,nth); ph_=np.linspace(0,2*np.pi,nph,endpoint=False)
    T,P=np.meshgrid(th,ph_,indexing='ij')
    d=np.stack([np.sin(T)*np.cos(P),np.sin(T)*np.sin(P),np.cos(T)],-1)
    Q=(c+(d*r)@np.linalg.inv(rec)).reshape(-1,3)
    ph.run_qpoints([list(q) for q in Q],with_eigenvectors=True)
    V=ph.qpoints.eigenvectors; nb=V.shape[1]
    V=V.reshape(nth,nph,nb,nb)[:,:,:,:NOCC]
    L=lambda A,B:(lambda m:m/np.abs(m))(np.linalg.det(np.einsum('...ji,...jk->...ik',A.conj(),B)))
    Vp=np.concatenate([V,V[:,:1]],axis=1)
    F=np.angle(L(Vp[:-1,:-1],Vp[:-1,1:])*L(Vp[:-1,1:],Vp[1:,1:])/(L(Vp[1:,:-1],Vp[1:,1:])*L(Vp[:-1,:-1],Vp[1:,:-1])))
    return F.sum()/(2*np.pi)
def wcc(ph,rec,c,r,nth=81,nph=120):
    ths=np.linspace(1e-3,np.pi-1e-3,nth); phis=np.linspace(0,2*np.pi,nph,endpoint=False)
    s=[]
    for t in ths:
        d=np.stack([np.sin(t)*np.cos(phis),np.sin(t)*np.sin(phis),np.cos(t)*np.ones(nph)],-1)
        Q=c+(d*r)@np.linalg.inv(rec)
        ph.run_qpoints([list(q) for q in Q],with_eigenvectors=True)
        V=ph.qpoints.eigenvectors[:,:,:NOCC]
        W=np.eye(NOCC,dtype=complex)
        for i in range(nph): W=W@(V[i].conj().T@V[(i+1)%nph])
        s.append(np.angle(np.linalg.det(W)))
    d=np.diff(np.array(s)); d=(d+np.pi)%(2*np.pi)-np.pi
    return d.sum()/(2*np.pi)

CASES=[('pristine',(0.0975234,0.0,0.1610282)),
       ('m1 (-1%)',(0.0637717,0.0,0.1294182)),
       ('p1 (+1%)',(0.2643381,0.0,0.3182116)),
       ('parent',  (0.1464927,0.0708493,0.0))]
SRC={'pristine':'pristine_mirror','m1 (-1%)':'m1_mirror','p1 (+1%)':'p1_mirror','parent':'parent_pristine'}
RAD=[0.0010,0.0020,0.0035,0.0050,0.0080,0.0120]
for lab,k in CASES:
    D='/home/user/HfO2-Orthorhombic-Polar-Project/source/'+SRC[lab]
    ph=get_ph(D); rec=np.linalg.inv(ph.primitive.cell).T*2*np.pi
    k=np.array(k)
    g=gap_at(ph,[list(k)],NOCC)[0]
    print('%-10s gap %.2e' % (lab,g), flush=True)
    print('   flux(60x120): ' + '  '.join('r=%.4f:%+.2f'%(r,flux(ph,rec,k,r)) for r in RAD), flush=True)
    print('   WCC (81x120): ' + '  '.join('r=%.4f:%+.2f'%(r,wcc(ph,rec,k,r)) for r in (0.0020,0.0035,0.0080)), flush=True)
