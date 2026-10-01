"""검증 (D): 노드 주변 2x2 유효 해밀토니안의 선형 전개.

D(k) 를 축퇴 부분공간에 사영해 H = f0*I + sum_a f_a sigma_a 로 쓰고,
f_a = sum_i v_ai dk_i 의 3x3 속도행렬 v 를 최소제곱으로 뽑는다.

  - chi = sign det(v)  -> Berry flux 와 완전히 다른 경로의 chirality
  - 선형성 잔차가 작아야 진짜 선형 교차(Weyl). 2차 접촉이면 v 가 특이해진다.
  - tilt = f0 의 기울기. |tilt| > |cone| 이면 type-II.
"""
import sys, numpy as np
sys.path.insert(0, __file__.rsplit('/',2)[0])
from weyl_scan import get_ph

SRC=sys.argv[1]; K0=np.array([float(x) for x in sys.argv[2:5]])
NB=17                      # 17,18 번 밴드 (0-based 16,17)
R=float(sys.argv[5]) if len(sys.argv)>5 else 1e-4
ph=get_ph(SRC)
dm=ph.dynamical_matrix

def D(k):
    dm.run(list(k))
    M=dm.dynamical_matrix
    return (M+M.conj().T)/2

D0=D(K0); w0,V0=np.linalg.eigh(D0)
V=V0[:,NB-1:NB+1]                       # 축퇴 2개
print('node k = (%.7f, %.7f, %.7f)'%tuple(K0))
print('band 17,18 eigenvalue  %.9f  %.9f   (차이 %.3e)'%(w0[NB-1],w0[NB],w0[NB]-w0[NB-1]))
FAC=ph.unit_conversion_factor
print('omega0 = %.6f THz'%(np.sqrt(abs(w0[NB-1]))*FAC))

sig=[np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]]),np.array([[1,0],[0,-1]],complex)]
# 유한차분 표본: 각 축 +-, 그리고 대각 방향 몇 개
dirs=[]
for i in range(3):
    e=np.zeros(3); e[i]=1; dirs += [e,-e]
for a,b in [(0,1),(0,2),(1,2)]:
    e=np.zeros(3); e[a]=1; e[b]=1; dirs += [e/np.sqrt(2), -e/np.sqrt(2)]
    e=np.zeros(3); e[a]=1; e[b]=-1; dirs += [e/np.sqrt(2), -e/np.sqrt(2)]
P0=V.conj().T @ D0 @ V; P0=(P0+P0.conj().T)/2
Y0=np.array([np.real(np.trace(P0))/2]+[np.real(np.trace(P0@s))/2 for s in sig])
X=[]; Y=[]
for d in dirs:
    dk=R*d
    P=V.conj().T @ D(K0+dk) @ V
    P=(P+P.conj().T)/2
    c0=np.real(np.trace(P))/2
    ca=[np.real(np.trace(P@s))/2 for s in sig]
    X.append(dk); Y.append(np.array([c0]+ca)-Y0)
X=np.array(X); Y=np.array(Y)
coef,res,rank,sv=np.linalg.lstsq(X,Y,rcond=None)     # (3 dk) x (4 comp)
tilt=coef[:,0]; v=coef[:,1:].T                        # v[a,i]
pred=X@coef
resid=np.abs(Y-pred).max()/np.abs(Y).max()
r2=1-((Y-pred)**2).sum()/((Y-Y.mean(0))**2).sum()
detv=np.linalg.det(v)
print('\n속도행렬 v (a=x,y,z sigma / i=k1,k2,k3), 단위 (THz^2)/rlu')
for a in range(3):
    print('   %+.6e  %+.6e  %+.6e'%tuple(v[a]))
print('det v = %+.6e   -> chi = %+d'%(detv, int(np.sign(detv))))
print('선형 피팅: 상대잔차 %.2e,  R^2 = %.8f'%(resid,r2))
print('v 의 특이값 %s'%np.array2string(np.linalg.svd(v,compute_uv=False),precision=4))
print('조건수 %.3e (특이하면 2차 접촉)'%np.linalg.cond(v))
print('\ntilt (f0 기울기) %s'%np.array2string(tilt,precision=4))
smin=np.linalg.svd(v,compute_uv=False).min()
print('|tilt| / (cone 최소축) = %.2f  -> %s'%(np.linalg.norm(tilt)/smin,
      'type-II' if np.linalg.norm(tilt)/smin>1 else 'type-I'))
