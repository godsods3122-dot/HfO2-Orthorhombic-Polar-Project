"""표면 측정 쪽: 노드 에너지에서 투영된 bulk 연속체가 표면 BZ 를 얼마나 덮는가.
SURFACE 카드 (1 0 0 / 0 1 0) -> 법선이 b3(편극축), 즉 k3 방향으로 투영."""
import sys, numpy as np
sys.path.insert(0,'scripts')
from weyl_scan import get_ph
THZ2MEV=4.135667696; E0=10.0868
ph=get_ph('source/parent_pristine',(2,2,2))
N1,N3=41,41
k1=np.linspace(-0.5,0.5,N1); k2=np.linspace(-0.5,0.5,N1); k3=np.linspace(-0.5,0.5,N3)
MIN=np.zeros((N1,N1,36)); MAX=np.zeros((N1,N1,36))
for i,a in enumerate(k1):
    qs=[[a,c,z] for c in k2 for z in k3]
    ph.run_qpoints(qs); f=ph.qpoints.frequencies.reshape(N1,N3,36)
    MIN[i]=f.min(axis=1); MAX[i]=f.max(axis=1)
    if i%10==0: print("  k1 %d/%d"%(i,N1),flush=True)
print("\n=== 표면 BZ 에서 노드 에너지 E0=%.4f THz 를 덮는 투영 bulk 밴드 ==="%E0)
cov=((MIN<=E0)&(MAX>=E0))
nb=cov.sum(axis=2)
print("  덮는 밴드 수 분포:")
for v in range(0,int(nb.max())+1):
    frac=(nb==v).mean()
    if frac>0: print("    %d개 밴드가 덮음 : 표면 BZ 의 %5.1f%%"%(v,100*frac))
print("  -> 진짜 gap(덮는 밴드 0개) 인 면적: %.1f%%"%(100*(nb==0).mean()))
print("  -> bulk 가 덮는 면적            : %.1f%%"%(100*(nb>0).mean()))
bands=sorted(set(np.where(cov.any(axis=(0,1)))[0]))
print("  관여 밴드: %s"%[b+1 for b in bands])
# 노드의 표면 투영 위치 주변
i=np.argmin(np.abs(k1-0.1464927)); j=np.argmin(np.abs(k2-0.0708493))
print("\n  노드 투영점 (k1,k2)=(%.3f,%.3f) 부근 5x5 의 덮는 밴드 수:"%(k1[i],k2[j]))
print(nb[i-2:i+3,j-2:j+3])
np.savez('/tmp/proj.npz',MIN=MIN,MAX=MAX,k1=k1,k2=k2)
