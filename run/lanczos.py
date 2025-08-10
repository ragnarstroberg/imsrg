from pyIMSRG import *
import numpy as np

cm=Commutator
gm=Generator()

## initialize the T and D^dagger

def htc(Haml, chi):
    ## generate a antihermit chi
    chi_d=chi*0
    chi_d = gm.GetVSEOM_ladder(chi, 1)
    chi_d.SetAntiHermitian()

    ht_plus= chi*0
    ht_minus= chi*0
    ## chi is hermitian, hplus is antihermitian

    ht_plus.SetAntiHermitian()
    ht_minus.SetHermitian()

    ht_plus = cm.Commutator(Haml, chi )
    ht_minus = cm.Commutator(Haml, chi_d )


    heom1= gm.GetVSEOM_ladder(ht_plus, 0)
    heom2= gm.GetVSEOM_ladder(ht_minus, 0)
   # hod=heom1
    hod = (heom1+heom2)/2

    return(hod)

def htc_vs(Haml, chi):


    chi_d=chi*0
    chi_d.SetAntiHermitian()
    chi_d = gm.GetVSEOM_ladder(chi, 1)

    ht_plus= chi*0
    ht_minus= chi*0
    ## chi is hermitian, hplus is antihermitian

    ht_plus.SetAntiHermitian()

    ht_minus.SetHermitian()

    ht_plus = cm.Commutator(Haml, chi )
    ht_minus = cm.Commutator(Haml, chi_d )

#    cm.comm111ss(Haml, chi  , ht_plus  )
#    cm.comm111ss(Haml, chi_d, ht_minus )
#    cm.comm121ss(Haml, chi  , ht_plus  )
#    cm.comm121ss(Haml, chi_d, ht_minus )
#    cm.comm222_pp_hhss(Haml, chi  , ht_plus  )
#    cm.comm222_pp_hhss(Haml, chi_d, ht_minus )
#    cm.comm222_phss(Haml, chi  , ht_plus  )
#    cm.comm222_phss(Haml, chi_d, ht_minus )
#    cm.comm221ss(Haml, chi  , ht_plus  )
#    cm.comm221ss(Haml, chi_d, ht_minus )
#    cm.comm122ss(Haml, chi  , ht_plus  )
#    cm.comm122ss(Haml, chi_d, ht_minus )
#    cm.comm222_pp_hh_221ss(Haml, chi  , ht_plus  )
#    cm.comm222_pp_hh_221ss(Haml, chi_d, ht_minus )


    heom1= gm.GetVSEOM_ladder(ht_plus, 0)
    heom2= gm.GetVSEOM_ladder(ht_minus, 0)

    hod=Haml*0
    hod.SetHermitian()
    hod = (heom1+heom2)/2
    #hod = (ht_plus+ht_minus)
    #hod=gm.GetVSEOM_ladder(hod,0)
    return(hod)


def Norm(T1, T2):
    return(gm.GetEOM_Overlap(T1,T2))


def Norm_vs_new(T1,T2):

    T1d=gm.GetVSEOM_ladder(T1,1)
    T2d=gm.GetVSEOM_ladder(T2,1)
    nop=T1*0
    nop1=T1*0
    nop2=T1*0
    nop3=T1*0
    nop4=T1*0
    nop.SetNonHermitian()
    nop1.SetAntiHermitian()
    nop2.SetHermitian()
    nop3.SetHermitian()
    nop4.SetAntiHermitian()

    nop1= cm.Commutator(T1,T2)
    nop2= cm.Commutator(T1d,T2)
    nop3= cm.Commutator(T1,T2d)
    nop4= cm.Commutator(T1d,T2d)

    nop=nop1-nop2+nop3-nop4
    nop=nop/4
    nm = gm.GetVSEOM_Overlap(nop)
    return(nm)

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

#    w = hv_func(haml,lanczos_vector[0])
#
#    ai=norm_func(w,lanczos_vector[0])
#    print('first step ai :', ai)

    

    for j in range(max_iter):

        w = hv_func(haml,lanczos_vector[j])

        ai=norm_func(w,lanczos_vector[j])
        print('ai at', j, ': ', ai )


        if(j>0):
            w=w-ai*lanczos_vector[j]-bj*lanczos_vector[j-1]
        else:
            w=w-ai*lanczos_vector[j]

        hall[j,j]=ai


        bj = np.sqrt(norm_func(w,w))

        if bj < 0.01 :
            print('bj is small')
            break
        lanczos_vector.append(w/bj)

        if(j<max_iter-1):
            hall[j,j+1]=bj
            hall[j+1,j]=bj
        print(j,ai,bj)
        if(j > 2 and j%1 == 0):
            e,v = np.linalg.eig(hall[0:j,0:j])
            e=np.sort(e)
            print("Energy on ", j , ' th iteration: ', e[0:state_want] )
            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.01):
                print('energy converge')
                break
            norm_e_old=norm_e_new

    print(hall)
    print('Lanczos converged with ', j, ' step')
    for k in range(state_want):
        print('E(',k,')= ', e[k])
    print('Lanczos converged with ', j, ' step')
    return(e[0:state_want], v[0:state_want,:],lanczos_vector)





def arnoldi_proc( hv_func, norm_func, haml, vi, max_iter,state_want):
    lanczos_vector = []
    hall = np.zeros([max_iter,max_iter])
    hall[0,0]=0.

    ## normalize it to 1
    nn=norm_func(vi,vi)
    print('initial norm vector', nn)
    vi=vi/np.sqrt(nn)
    lanczos_vector.append(vi)
    norm_e_old=-10000
    norm_e_new=-1000
    e=np.zeros(state_want)

    for j in range(max_iter-1):
        w = hv_func(haml,lanczos_vector[j])
        ai=norm_func(w,lanczos_vector[j])
        hall[j,j]=ai
        w=w-ai*lanczos_vector[j]

        if(j>0):
            for i in range(j):
                bj=norm_func(w,lanczos_vector[i])
                hall[j,i]=bj
                w=w-bj*lanczos_vector[i]

        ## generate the new vector
        bj = norm_func(w,w)
        if bj < 0.1:
            print('bj is small: ', bj, j)
            break
        print(j,'th step: ',bj)
        w=w/np.sqrt(bj)
        lanczos_vector.append(w)

        for m in range(j):
            vl=lanczos_vector[j]
            vm=lanczos_vector[m]
            vec=hv_func(haml, vm)
            nmm=norm_func(vl,vec)
            hall[m,j] =nmm
        if(j >= state_want ):
            #hsub=(hsub+hsub.T)/2
            e,v = np.linalg.eig(hall[0:j+1,0:j+1])
            e=np.sort(e)
            print("Energy on ", j , ' th iteration: ', e[0:state_want] )
            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.01):
                print('energy converge')
                break
            norm_e_old=norm_e_new

    print(hall)
    print('Arnoldi converged with ', j, ' step')
    #for k in range(state_want):
    #    print('E(',k,')= ', e[k])
    #print('Arnoldi converged with ', j, ' step')

    return(e, v,lanczos_vector)
    #return(lanczos_vector)
