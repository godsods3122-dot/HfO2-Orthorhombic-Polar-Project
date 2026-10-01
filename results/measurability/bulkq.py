"""bulk 가 노드 측정을 왜 어렵게 하는가 — 정량 분석."""
import sys, numpy as np
sys.path.insert(0,'scripts')
from weyl_scan import get_ph
THZ2MEV=4.135667696
E0=10.0868                      # 노드 주파수 THz
ph=get_ph('source/parent_pristine',(2,2,2))
B=np.linalg.inv(ph.primitive.cell).T*2*np.pi; Binv=np.linalg.inv(B)
NODE=np.array([0.1464927,0.0708493,0.0]); NC=NODE@B
print("노드 E0 = %.4f THz = %.3f meV,  cartesian k = %s"%(E0,E0*THZ2MEV,np.round(NC,5)))

# ---------- 1. 상태밀도: 이 에너지가 얼마나 붐비는가 ----------
z=np.load('/tmp/bg.npz'); F=z['F']; nq=F.shape[0]
V=abs(np.linalg.det(B))         # BZ 부피 (1/Å^3)
print("\n=== 1. 에너지 붐빔 (BZ 표본 %d q × 36모드, BZ 부피 %.4f 1/Å³) ==="%(nq,V))
for dE_meV in (0.5,1.0,1.5,3.0):
    dE=dE_meV/THZ2MEV
    sel=np.abs(F-E0)<dE
    # 단위포당 상태밀도 g(E0) [states/meV/cell]
    g=sel.sum()/nq/(2*dE_meV)
    print("  ±%.1f meV: 모드 %6d개,  q점의 %5.1f%%,  g = %.4f states/meV/cell,  관여밴드 %s"
          %(dE_meV,sel.sum(),100*sel.any(axis=1).mean(),g,
            [int(b)+1 for b in sorted(set(np.where(sel)[1]))]))
# 전 스펙트럼 대비
hist,edges=np.histogram(F.ravel(),bins=120,range=(0,23))
gmax=hist.max()/nq/((edges[1]-edges[0])*THZ2MEV)
i0=np.searchsorted(edges,E0)-1
g0=hist[i0]/nq/((edges[1]-edges[0])*THZ2MEV)
print("  g(E0) = %.4f,  스펙트럼 최대 g = %.4f  ->  E0 는 최대치의 %.0f%%"%(g0,gmax,100*g0/gmax))

# ---------- 2. 노드는 점인가 선인가 ----------
print("\n=== 2. 분해능 안에서 노드의 실제 모양 ===")
def gap(kc):
    ph.run_qpoints([list(kc@Binv)]); f=ph.qpoints.frequencies[0]; return float(f[17]-f[16])
dirs={'cart x (면내)':[1,0,0],'cart y (편극축, 느린 축)':[0,1,0],'cart z (면내)':[0,0,1]}
print("  방향별로 gap 이 분해능에 닿는 거리")
for lbl,d in dirs.items():
    d=np.array(d,float)
    row=[]
    for target_meV in (0.1,1.0,3.0):
        tgt=target_meV/THZ2MEV
        lo,hi=1e-5,2.0
        for _ in range(60):
            mid=(lo+hi)/2
            if gap(NC+mid*d)<tgt: lo=mid
            else: hi=mid
        row.append(lo)
    print("   %-26s  0.1meV: %8.4f    1meV: %8.4f    3meV: %8.4f  (1/Å)"%(lbl,*row))
