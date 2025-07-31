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
    #hod=heom1
    hod = (heom1+heom2)/2

    return(hod)

def htc_vs(Haml, chi):


    ## generate a antihermit chi
    chi_d = gm.GetVSEOM_ladder(chi, 1)

    ht_plus= chi*0
    ht_minus= chi*0
    ## chi is hermitian, hplus is antihermitian

    ht_plus.SetAntiHermitian()

    ht_minus.SetHermitian()

    ht_plus = cm.Commutator(Haml, chi )

    ht_minus = cm.Commutator(Haml, chi_d )


    heom1= gm.GetVSEOM_ladder(ht_plus, 0)
    heom2= gm.GetVSEOM_ladder(ht_minus, 0)
    hod=heom2
    hod = (heom1+heom2)/2

    return(hod)


def Norm(T1, T2):
    return(gm.GetEOM_Overlap(T1,T2))

def Norm_vs(T1,T2):

    T0=gm.GetVSEOM_ladder(T2,1)
    nop=T1*0
    nop.SetHermitian()
    nop=cm.Commutator(T1,T0)
    nop=nop/2

    nm = gm.GetVSEOM_Overlap(nop)
    return(nm)

import numpy as np

def lanczos_proc( hv_func, norm_func, haml, vi, max_iter,state_want):
    lanczos_vector = []
    hall = np.zeros([max_iter,max_iter])
    hall[0,0]=0.

    ## normalize it to 1
    nn=norm_func(vi,vi)
    print('initial norm vector', nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)
    norm_e_old=-1000
    norm_e_new=-1000

    w = hv_func(haml,lanczos_vector[0])

    ai=norm_func(w,lanczos_vector[0])
    print('first step ai :', ai)

    

    for j in range(max_iter):

        w = hv_func(haml,lanczos_vector[j])

        ai=norm_func(w,lanczos_vector[j])


        if(j>0):
            w=w-ai*lanczos_vector[j]-bj*lanczos_vector[j-1]
        else:
            w=w-ai*lanczos_vector[j]

        hall[j,j]=ai


        bj = np.sqrt(norm_func(w,w))

        if bj < 0.0000001 :
            break
        lanczos_vector.append(w/bj)

        if(j<max_iter-1):
            hall[j,j+1]=bj
            hall[j+1,j]=bj
        print(j,ai,bj)
        if(j > 2 and j%2 == 0):
            e,v = np.linalg.eig(hall[0:j,0:j])
            e=np.sort(e)
            print("Energy on ", j , ' th iteration: ', e[0:state_want] )
            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.01):
                break
            norm_e_old=norm_e_new

    print('Lanczos converged with ', j, ' step')
    for k in range(state_want):
        print('E(',k,')= ', e[k])
    print('Lanczos converged with ', j, ' step')
    return(e[0:state_want], v[0:state_want,:],lanczos_vector)



