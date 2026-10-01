#!/usr/bin/env python3
"""
mode_dipole_and_eps.py — 극성(=ferron) 성분의 크기를 정량화한다.

(1) 이온(=ferron) 이 유전상수에 기여하는 양, 모드별로.
    eps^ion_ab = (1/(eps0*V)) sum_m  Zt_m,a Zt_m,b / w_m^2 ,  Zt = sum_j e Z*_j e_j/sqrt(m_j)
(2) BZ 전체에서 모드 쌍극자 |p| 분포 — 노드가 극성에서 특별한가?"""
import sys, numpy as np
sys.path.insert(0,'scripts')
from weyl_scan import get_ph
E=1.602176634e-19; EPS0=8.8541878128e-12; U=1.66053906660e-27; HBAR=1.054571817e-34
ph=get_ph('source/parent_pristine',(2,2,2))
prim=ph.primitive; m=np.array(prim.masses); Z=np.array(ph.nac_params['born'])
V=abs(np.linalg.det(prim.cell))*1e-30
B=np.linalg.inv(prim.cell).T*2*np.pi
ph.run_qpoints([[0.,0.,0.]],with_eigenvectors=True)
f=ph.qpoints.frequencies[0]; ev=ph.qpoints.eigenvectors[0]
eps=np.zeros((3,3)); contrib=[]
for mo in range(36):
    w=f[mo]*1e12*2*np.pi
    if f[mo]<0.5: continue
    e=ev[:,mo].reshape(-1,3)
    Zt=E*np.einsum('jab,jb->a',Z,e.real/np.sqrt(m*U)[:,None])
    d=np.outer(Zt,Zt)/(EPS0*V*w**2)
    eps+=d; contrib.append((mo,f[mo],np.trace(d)/3,d))
print("=== 이온(ferron) 유전 기여  eps^ion  (eps_inf 는 별도) ===")
print(np.round(eps,3))
print("  대각: (%.2f, %.2f, %.2f)   평균 %.2f"%(*np.diag(eps),np.trace(eps)/3))
tot=np.trace(eps)/3
print("\n  모드별 기여 (평균 대각 기준, 상위 10개)")
for mo,fr,c,_ in sorted(contrib,key=lambda t:-t[2])[:10]:
    print("   band %2d  %7.3f THz   %7.3f  (%5.1f%%)"%(mo+1,fr,c,100*c/tot))
near=[c for mo,fr,c,_ in contrib if abs(fr-10.0868)<0.5]
print("\n  노드 에너지 10.09±0.5 THz 모드들의 합계 기여: %.3f (%.1f%%)"%(sum(near),100*sum(near)/tot))

print("\n=== BZ 전체 모드 쌍극자 |p| 분포 (e·Å) ===")
N=12; g=(np.arange(N)+0.5)/N-0.5
QS=np.array([[x,y,z] for x in g for y in g for z in g])
P=np.zeros((len(QS),36)); C=np.zeros((len(QS),36)); F=np.zeros((len(QS),36))
for s in range(0,len(QS),300):
    q=QS[s:s+300]
    ph.run_qpoints([list(x) for x in q],with_eigenvectors=True)
    F[s:s+300]=ph.qpoints.frequencies
    e=ph.qpoints.eigenvectors.transpose(0,2,1).reshape(len(q),36,12,3)
    Zv=np.einsum('jab,qmjb->qma',Z,e/np.sqrt(m)[None,None,:,None])
    w=np.maximum(ph.qpoints.frequencies,1e-3)
    A=np.sqrt(HBAR/(2*2*np.pi*w*1e12*U))*1e10
    p=A[:,:,None]*Zv
    P[s:s+300]=np.linalg.norm(p,axis=2)
    a,b=p.real,p.imag
    cr=np.linalg.norm(np.cross(a,b),axis=2); den=(a*a).sum(2)+(b*b).sum(2)
    C[s:s+300]=np.where(den>1e-14,2*cr/np.maximum(den,1e-30),0)
def st(x,l): print("  %-26s 중앙값 %.4f 평균 %.4f 90%%분위 %.4f 최대 %.4f"%(l,np.median(x),x.mean(),np.percentile(x,90),x.max()))
ok=F>0.5
st(P[ok],"|p| 전 밴드 x 전 BZ")
st(P[:,16:18].ravel(),"|p| band 17,18")
sel=(np.abs(F-10.0868)<0.10)&ok
st(P[sel],"|p| 노드 에너지 ±0.10")
print()
st(C[ok],"쌍극자 원형도 전체")
st(C[:,16:18].ravel(),"쌍극자 원형도 b17,18")
st(C[sel],"쌍극자 원형도 노드에너지")
