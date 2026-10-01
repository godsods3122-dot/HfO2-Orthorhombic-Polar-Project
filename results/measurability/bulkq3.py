import sys, numpy as np
sys.path.insert(0,'scripts')
from weyl_scan import get_ph
THZ2MEV=4.135667696; E0=10.0868
ph=get_ph('source/parent_pristine',(2,2,2))
B=np.linalg.inv(ph.primitive.cell).T*2*np.pi; Binv=np.linalg.inv(B)
NC=np.array([0.1464927,0.0708493,0.0])@B
def bands(kc):
    ph.run_qpoints([list(kc@Binv)]); return ph.qpoints.frequencies[0]
print("=== 2. gap 이 분해능에 닿으려면 노드에서 얼마나 멀어져야 하나 (스캔, 첫 도달) ===")
res={}
for lbl,d in (('cart x (면내)',[1,0,0]),('cart y (편극축)',[0,1,0]),('cart z (면내)',[0,0,1])):
    d=np.array(d,float); ts=np.linspace(1e-4,0.60,240)
    g=np.array([bands(NC+t*d)[17]-bands(NC+t*d)[16] for t in ts])*THZ2MEV
    row=[]
    for tgt in (0.1,0.5,1.5,3.0):
        i=np.argmax(g>=tgt)
        row.append(ts[i] if g.max()>=tgt else np.nan)
    res[lbl]=row
    print("  %-16s 0.1meV %7.4f | 0.5meV %7.4f | 1.5meV %7.4f | 3.0meV %7.4f   (1/Å)"%(lbl,*row))
print("\n=== 3. 그 흐릿한 영역 안에 bulk 가 몇 개나 더 있나 ===")
print("  ΔE     R_blur(x,y,z)         BZ 부피비     그 안에서 ±ΔE 에 드는 밴드")
V=abs(np.linalg.det(B))
for k,tgt in enumerate((0.5,1.5,3.0)):
    Rx,Ry,Rz=res['cart x (면내)'][k+1],res['cart y (편극축)'][k+1],res['cart z (면내)'][k+1]
    if np.isnan(Rx) or np.isnan(Ry) or np.isnan(Rz): print("  %.1f meV  (BZ 안에서 도달 못함)"%tgt); continue
    vol=4/3*np.pi*Rx*Ry*Rz
    # 그 타원체 안을 성기게 샘플해 어떤 밴드가 E0±ΔE 안에 있는지
    rng=np.random.default_rng(0); n=600
    u=rng.normal(size=(n,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
    r=rng.random(n)**(1/3)
    P=NC+u*r[:,None]*np.array([Rx,Ry,Rz])
    ph.run_qpoints([list(p@Binv) for p in P]); F=ph.qpoints.frequencies
    sel=np.abs(F-E0)<tgt/THZ2MEV
    bl=[int(b)+1 for b in sorted(set(np.where(sel)[1]))]
    print("  %.1f meV  (%.4f, %.4f, %.4f)   %6.2f%%      %s  (평균 %.2f개/q점)"
          %(tgt,Rx,Ry,Rz,100*vol/V,bl,sel.sum()/n))
