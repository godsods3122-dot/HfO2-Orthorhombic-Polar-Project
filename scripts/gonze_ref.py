"""phonopy의 Gonze dipole-dipole을 순수 numpy로 재구현 (Fortran 이식 전 검증용)"""
import numpy as np, itertools

def build_G_list(rec_lat, G_cutoff, g_rad=100):
    # rec_lat: columns are b1,b2,b3 (phonopy: self._rec_lat = pcell.reciprocal_lattice)
    for _g in range(g_rad, 0, -1):
        found=False
        for a,b,c in itertools.product((-1,0,1),repeat=3):
            if (a,b,c)==(0,0,0): continue
            if np.linalg.norm(rec_lat @ [a,b,c])*_g < G_cutoff:
                g_rad=_g+1; found=True; break
        if found: break
    rng=np.arange(-g_rad,g_rad+1)
    idx=np.array(list(itertools.product(rng,rng,rng)))
    G=idx @ rec_lat.T
    n2=(G**2).sum(axis=1)
    return G[n2 < G_cutoff**2]

def recip_dd(q_cart, q_dir_cart, G_list, born, eps, pos_cart, Lambda, tol=1e-5):
    """pos_cart: phonopy의 self._pcell.scaled_positions? -> C코드의 pos는 scaled positions"""
    n=len(born); L2=4*Lambda**2
    KK=np.zeros((len(G_list),3,3))
    for g,G in enumerate(G_list):
        qK=G+q_cart; norm=np.sqrt((qK**2).sum())
        if norm < tol:
            if q_dir_cart is None: continue
            dp=q_dir_cart @ eps @ q_dir_cart
            KK[g]=np.outer(q_dir_cart,q_dir_cart)/dp
        else:
            dp=qK @ eps @ qK
            KK[g]=np.outer(qK,qK)/dp*np.exp(-dp/L2)
    dd_part=np.zeros((n,3,n,3),dtype=complex)
    for g,G in enumerate(G_list):
        ph=np.exp(2j*np.pi*(pos_cart @ G))       # phase = 2pi*(pos_i-pos_j).G
        M=np.outer(ph,ph.conj())                  # [i,j] = exp(2pi i (pos_i-pos_j).G)
        dd_part += KK[g][None,:,None,:]*M[:,None,:,None]
    # multiply Born charges: dd[i,a,j,b] = sum_{k,l} Z[i][k][a] dd_part[i,k,j,l] Z[j][l][b]
    dd=np.einsum('ika,ikjl,jlb->iajb', born, dd_part, born)
    return dd
