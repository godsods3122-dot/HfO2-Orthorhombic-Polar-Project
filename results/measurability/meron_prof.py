"""메론의 위치·크기·특성을 극좌표 프로파일로 잰다.

메론은 (polarity p, vorticity w, helicity gamma) 로 특징지어지고
    Q = (w/2) [ n_p(0) - n_p(inf) ]
이므로, 반지름별 <n_polar> 와 면내 감김수만 재면 Q 까지 나온다.
격자 전체를 도는 것보다 훨씬 싸다.
"""
import sys, numpy as np
sys.path.insert(0,'scripts'); sys.path.insert(0,'scripts/figs')
import matplotlib; matplotlib.use('Agg')
from fig9_berry_field import Berry
b=Berry('source/parent_pristine'); nb=b.n
NODE=np.array([0.1464927,0.0708493,0.0])*nb
NANG=48
ang=np.linspace(0,2*np.pi,NANG,endpoint=False)

def ring(r,delta,h):
    P=np.zeros((NANG,3))
    P[:,0]=NODE[0]+r*np.cos(ang); P[:,1]=NODE[1]+r*np.sin(ang); P[:,2]=delta
    O=b.omega([0.,0.,0.],P,h,comps=(0,1,2))
    nn=O/np.linalg.norm(O,axis=1,keepdims=True)
    npol=nn[:,2]
    # 면내 성분의 감김수 (vorticity) 와 helicity
    ph=np.arctan2(nn[:,1],nn[:,0])
    d=np.diff(np.concatenate([ph,ph[:1]])); d=(d+np.pi)%(2*np.pi)-np.pi
    w=d.sum()/(2*np.pi)
    # helicity: 면내 방향 - 방위각.  0 = 방사(바깥), pi = 방사(안쪽), +-pi/2 = 소용돌이
    g=np.angle(np.mean(np.exp(1j*(ph-ang))))
    return npol.mean(), npol.std(), w, g

print("노드 면내 위치 (환산) = (0.1464927, 0.0708493),  즉 (k1,k2) 그림자")
for delta in (0.0004, 0.002):
    print("\n===== 평면 k_polar = %+.4f  (스케일단위 1/Å) ====="%delta)
    print("   r(1/Å)    <n_polar>   흩어짐   감김수 w   helicity(도)")
    for r in (1e-4,2e-4,4e-4,8e-4,1.6e-3,3.2e-3,6.4e-3,1.28e-2,2.56e-2):
        m,s,w,g=ring(r,delta,h=min(r/4,2e-4))
        print("  %.5f   %+.4f    %.4f    %+.3f      %+7.1f"%(r,m,s,w,np.degrees(g)))
