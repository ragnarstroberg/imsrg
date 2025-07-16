from pyIMSRG import *
import numpy as np

cm=Commutator
gm=Generator()

## initialize the T and D^dagger

def htc(Haml, chi):

    ## generate a antihermit chi
    chi_d = gm.GetEOM_ladder(chi, 1)

    ht_plus= chi*0
    ht_minus= chi*0
    ## chi is hermitian, hplus is antihermitian

    ht_plus.SetAntiHermitian()
    ht_minus.SetHermitian()

    ht_plus = cm.Commutator(Haml, chi )

    ht_minus = cm.Commutator(Haml, chi_d )


    heom1= gm.GetEOM_ladder(ht_plus, 0)

    heom2= gm.GetEOM_ladder(ht_minus, 0)

    hod = (heom1+heom2)/2

    hod.SetHermitian()

    return(hod)




def Norm(T1, T2):
    return(gm.GetEOM_Overlap(T1,T2))

import numpy as np

def lanczos_proc( hv_func, norm_func, haml, vi, ndim):
    lanczos_vector = []
    hall = np.zeros([ndim,ndim])
    hall[0,0]=0.

    ## normalize it to 1
    nn=norm_func(vi,vi)
    print(nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)

    for j in range(ndim):

        w = hv_func(haml,lanczos_vector[j])

        ai=norm_func(w,lanczos_vector[j])


        if(j>0):
            w=w-ai*lanczos_vector[j]-bj*lanczos_vector[j-1]
        else:
            w=w-ai*lanczos_vector[j]

        hall[j,j]=ai


        bj = np.sqrt(norm_func(w,w))
        if bj < 0.00001 :
            break
        lanczos_vector.append(w/bj)

        if(j<ndim-1):
            hall[j,j+1]=bj
            hall[j+1,j]=bj
        #print(j,ai,bj)
    print(hall)
    e,v = np.linalg.eig(hall[0:j,0:j])

    return(e,v,lanczos_vector)

