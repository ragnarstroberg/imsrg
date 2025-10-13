
def Norm_vs_new(T1,T2,rdm):
    T1d=T1*0.
    T1d.SetAntiHermitian()
    T2d=T1*0.
    T2d.SetAntiHermitian()

    T1d=gm.GetVSEOM_ladder(T1,1)
    T2d=gm.GetVSEOM_ladder(T2,1)

    rst=0.

    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1,T2)
    rst += gm.GetVSEOM_Overlap(nop)


    nop=T1*0
    nop.SetHermitian()
    nop= cm.Commutator(T1,T2d)
    rst += gm.GetVSEOM_Overlap(nop)


    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1d,T2)
    rst -= gm.GetVSEOM_Overlap(nop)

    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1d,T2d)
    rst -= gm.GetVSEOM_Overlap(nop)
    return(rst/4.)
