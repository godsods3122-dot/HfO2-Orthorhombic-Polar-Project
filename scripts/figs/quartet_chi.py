"""4점 궤도 전부의 chirality 를 직접 계산 — 대칭 논증에 의존하지 않는다.

Pca2_1 의 k 공간 작용: C2(편극축, proper -> χ 보존), m_a 와 m_b (improper -> χ 반전).
따라서 1·3 사분면이 한 부호, 2·4 가 반대 부호여야 한다. 그걸 실측으로 확인한다.
"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/', 2)[0])
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from axes import roles, to_native
from weyl_scan import get_ph, gap_at
NOCC=17
def flux(ph,rec,c,r,nth=48,nph=96):
    th=np.linspace(0,np.pi,nth); pp=np.linspace(0,2*np.pi,nph,endpoint=False)
    T,P=np.meshgrid(th,pp,indexing='ij')
    d=np.stack([np.sin(T)*np.cos(P),np.sin(T)*np.sin(P),np.cos(T)],-1)
    Q=(c+(d*r)@np.linalg.inv(rec)).reshape(-1,3)
    ph.run_qpoints([list(q) for q in Q],with_eigenvectors=True)
    V=ph.qpoints.eigenvectors; nb=V.shape[1]
    V=V.reshape(nth,nph,nb,nb)[:,:,:,:NOCC]
    L=lambda A,B:(lambda m:m/np.abs(m))(np.linalg.det(np.einsum('...ji,...jk->...ik',A.conj(),B)))
    Vp=np.concatenate([V,V[:,:1]],axis=1)
    F=np.angle(L(Vp[:-1,:-1],Vp[:-1,1:])*L(Vp[:-1,1:],Vp[1:,1:])/(L(Vp[1:,:-1],Vp[1:,1:])*L(Vp[:-1,:-1],Vp[1:,:-1])))
    return F.sum()/(2*np.pi)

CASES=[('-0.8 %  m1_mirror',      'm1_mirror',      (0.1294182,0.0637717)),
       ('unstrained pristine',    'pristine_mirror',(0.1610282,0.0975234)),
       ('+1 %   p1_mirror',       'p1_mirror',      (0.3182116,0.2643381))]
print('4점 궤도 전부 직접 계산 (표준 k_a,k_b,  k_c=0).  자체 Berry flux 규약.')
for lab,src,(A,B) in CASES:
    D = __file__.rsplit('/', 3)[0] + '/source/' + src
    r2n,pol,mir=roles(D); ph=get_ph(D)
    rec=np.linalg.inv(ph.primitive.cell).T*2*np.pi
    print('\n%s'%lab)
    for sx in (1,-1):
        for sy in (1,-1):
            q=np.array(to_native((sx*A, sy*B, 0.0), r2n))
            g=gap_at(ph,[list(q)],NOCC)[0]
            chis=[flux(ph,rec,q,r) for r in (0.002,0.005)]
            print('   (%+.5f, %+.5f)  gap %.2e   chi = %s'
                  %(sx*A,sy*B,g,' / '.join('%+.2f'%c for c in chis)), flush=True)
