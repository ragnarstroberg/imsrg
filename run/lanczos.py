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

def Norm_vs_new2(T1,T2):
    T1d=T1*0.
    T1d.SetAntiHermitian()
    T2d=T1*0.
    T2d.SetAntiHermitian()
    nop2=T1*0.
    nop2.SetNonHermitian()

    T1d=gm.GetVSEOM_ladder(T1,1)
    T2d=gm.GetVSEOM_ladder(T2,1)

    rst=0.

    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1,T2)
    nop2+=nop


    nop=T1*0
    nop.SetHermitian()
    nop= cm.Commutator(T1,T2d)
    nop2+=nop


    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1d,T2)
    nop2-=nop

    nop=T1*0
    nop.SetNonHermitian()
    nop= cm.Commutator(T1d,T2d)
    nop2-=nop

    rst = gm.GetVSEOM_Overlap(nop2)
    return(rst/4.)


def Norm_vs_new3(T1,T2):
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

def Norm_vs_new(T1,T2):
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

def Norm3(t1,t2,haml,ms):

    rank_j, parity, rank_Tz, particle_rank = 0,0,0,3
    t3 = Operator(ms,rank_j, parity, rank_Tz, particle_rank)
    t3.ThreeBody.SetMode("pn")
    t1d=t1*0.
    t2d=t1*0.
    t1d.SetAntiHermitian()
    t2d.SetAntiHermitian()
    
    t1d=gm.GetVSEOM_ladder(t1,1)
    t2d=gm.GetVSEOM_ladder(t2,1)
    rst=0.

    nop=t1*0.0
    t3*=0.
    nop.SetHermitian()
    t3.SetAntiHermitian()
    cm.comm223ss(haml,t2,t3)
    cm.comm232ss(t1,t3,nop)
    cm.comm231ss(t1,t3,nop)
    cm.comm132ss(t1,t3,nop)
    rst=gm.GetVSEOM_Overlap(nop)

    nop=t1*0.0
    t3*=0.
    nop.SetAntiHermitian()
    t3.SetHermitian()
    cm.comm223ss(haml,t2d,t3)
    cm.comm232ss(t1,t3,nop)
    cm.comm231ss(t1,t3,nop)
    cm.comm132ss(t1,t3,nop)
    rst+=gm.GetVSEOM_Overlap(nop)
    

    nop=t1*0.0
    t3*=0.
    nop.SetAntiHermitian()
    t3.SetAntiHermitian()
    cm.comm223ss(haml,t2,t3)
    cm.comm232ss(t1d,t3,nop)
    cm.comm231ss(t1d,t3,nop)
    cm.comm132ss(t1d,t3,nop)
    rst-=gm.GetVSEOM_Overlap(nop)


    nop=t1*0.0
    t3*=0.
    nop.SetHermitian()
    t3.SetHermitian()
    cm.comm223ss(haml,t2d,t3)
    cm.comm232ss(t1d,t3,nop)
    cm.comm231ss(t1d,t3,nop)
    cm.comm132ss(t1d,t3,nop)
    rst-=gm.GetVSEOM_Overlap(nop)
    return(rst/4)



import numpy as np

def lanczos_proc( hv_func, norm_func, haml, vi, max_iter,state_want, ms):
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
        #ai-=Norm3(lanczos_vector[j],lanczos_vector[j],haml,ms)
        print('ai at', j, ': ', ai )


        if(j>0):
            w=w-ai*lanczos_vector[j]-bj*lanczos_vector[j-1]
        else:
            w=w-ai*lanczos_vector[j]

        hall[j,j]=ai

        bj = np.sqrt(norm_func(w,w))


        w1=w/bj
        w2=hv_func(haml,w1)
        nm2=norm_func(lanczos_vector[j],w2)
        nm3=Norm3(lanczos_vector[j],w1,haml, ms)
        bjj=nm2+nm3
        print('Beta, bj: ', bjj,bj)
        

        #bj = np.sqrt(norm_func(w,w))

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





def arnoldi_proc( hv_func, norm_func, haml, vi, max_iter,state_want,ms):
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
        dai=Norm3(lanczos_vector[j],lanczos_vector[j],haml, ms)
        ai=norm_func(w,lanczos_vector[j])
        hall[j,j]=ai+dai
        w=w-ai*lanczos_vector[j]

        if(j>0):
            for i in range(j):
                bj=norm_func(w,lanczos_vector[i])
#                hall[j,i]=bj
                w=w-bj*lanczos_vector[i]

        ## generate the new vector
        bj = norm_func(w,w)
        if bj < 0.1:
            print('bj is small: ', bj, j)
            break
        w=w/np.sqrt(bj)
        lanczos_vector.append(w)
        nmm=100
        nm3=100
        for m in range(j):
            vl=lanczos_vector[j]
            vm=lanczos_vector[m]
            vec=vl*0.0
            vec.SetHermitian()

            vec=hv_func(haml, vm)
            nmm=norm_func(vl,vec)
            nm3=Norm3(vl,vm,haml, ms)
            hall[m,j] =nmm+nm3
            hall[j,m] =nmm+nm3
        if nmm+nm3 < 0.1:
            print('coupling bj is small: ', bj, j)
            break
        print(j,'th step: ',bj,np.sqrt(bj),nmm+nm3, j+1,len(lanczos_vector))

        if(j >= state_want ):
            #hsub=(hall+hall.T)/2
            hsub=hall
          #  print(hsub[0:j+1,0:j+1])
            e,v = np.linalg.eig(hsub[0:j+1,0:j+1])
            idx = np.argsort(e)   
            e = e[idx]
            vs = v[:,idx]
            vec_a=lanczos_vector[j]*0.
            for kk in range(len(e)):
                vec_a = vec_a + lanczos_vector[kk]*vs[kk,0]
            vec_b=hv_func(haml, vec_a)
            nmm=norm_func(vec_a,vec_b)
            nm3=Norm3(vec_a,vec_a,haml, ms)
            nm=norm_func(vec_a,vec_a)
            print('energy check: ', nmm+nm3,nm)


            print("Energy on ", j , ' th iteration: ', e[0:state_want] )
            norm_e_new=0.0
            for k in range(state_want):
                norm_e_new += e[k]*e[k]
            if(abs(norm_e_new-norm_e_old) < 0.01):
                print('energy converge')
                break
            norm_e_old=norm_e_new
    hsub=hall[0:j,0:j]
    e,v = np.linalg.eig(hsub[0:j,0:j])
    e=np.sort(e)

    print('Arnoldi converged with ', j, ' step')
    #for k in range(state_want):
    #    print('E(',k,')= ', e[k])
    #print('Arnoldi converged with ', j, ' step')

    return(e, v,lanczos_vector)
    #return(lanczos_vector)
