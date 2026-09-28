"""특이점을 피해 링으로: 노드를 지나는 평면에서 면수직 Omega 가 정말 0 인가."""
import sys, numpy as np
sys.path.insert(0,'scripts'); sys.path.insert(0,'scripts/figs')
import matplotlib; matplotlib.use('Agg')
from fig9_berry_field import Berry
b=Berry('source/parent_pristine'); nb=b.n
R=np.array(b.B)/nb[:,None]
NODE=np.array([0.1464927,0.0708493,0.0])@np.array(b.B)
NA=48; ang=np.linspace(0,2*np.pi,NA,endpoint=False)
print("  평면       r(1/Å)   <|Omega_면내|>    <|Omega_y|>     비")
for delta in (0.0, 0.0004, 0.01):
    for r in (0.002, 0.01, 0.04):
        P=np.zeros((NA,3))
        P[:,0]=NODE[0]+r*np.cos(ang); P[:,1]=NODE[1]+delta; P[:,2]=NODE[2]+r*np.sin(ang)
        O=b.omega([0.,0.,0.],P@R.T,min(r/4,2e-3),comps=(0,1,2))@R
        ip=np.linalg.norm(O[:,[0,2]],axis=1).mean(); oo=np.abs(O[:,1]).mean()
        print("  ky%+0.4f  %.3f   %12.4e   %12.4e   %.3e"%(delta,r,ip,oo,oo/ip))
