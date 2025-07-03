import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import copy
from numpy import inf
from SIVZ_functions import *
from scipy.signal import find_peaks
from generic_functions import *

if __name__ == '__main__':
    indice=int(sys.argv[1]) # 0= small diatom, 1= picoeukaryoe, 2= symechococcus, 3=prochlorococcus
    eff=sys.argv[2] # loss of infectivity (0)
    otype=sys.argv[3] # ocean type: 0, oligotrophic, mesotrophic or upwelling
    dz2=float(sys.argv[4]) # qudratic mortality of Z
    m2=float(sys.argv[5]) # quadratic mortality of V
    alph=float(sys.argv[6]) # strenght of grazing on I
    pred_type=str(sys.argv[7]) # 0 (=DARWIN parameterization of Z) or encounter_model= biophysical encounter rate
    big=str(sys.argv[8]) # 0 =smalleest cells, 1= biggest cells

    # print arguments
    print(sys.argv)

    # load life history traits
    Vols=load_vector('../trait_data/Vs_5.txt', sep=' ')
    Ncs=load_vector('../trait_data/Nc_dutkiewicz_5.txt', sep=' ')
    betas=load_vector('../trait_data/model_burst_size_nn-gam.txt', sep=' ')
    lps=load_vector('../trait_data/model_latent_period_nn-gam.txt', sep=' ')
    mu_max=load_vector('../trait_data/mumax_dutkiewicz_5.txt', sep=' ')

    if indice in [2,3]:
        rz=2.5
    elif indice in [0,1]:
        rz=5
    # nitrogen quota of the grazer
    Qz=Q_grazer(rz)

    # nitrogen quota of the virus
    r_virus=[20,80,35,35]
    Qvs=[Q_virus(r) for r in r_virus]

    Qps=[]
    for i in range(100):
        Qps.append(Q_diatom(Vols[i]))
    for i in range(100,200):
        Qps.append(Q_eukaryotes(Vols[i]))
    for i in range(200,400):
        Qps.append(Q_cyanobacteria(Vols[i]))


    N=4 # number of algae

    if N==4:
        # choose life history trait of either smallest or biggest cells
        if big=='0':
            ind1=np.argmin(np.array(Qps[0:100]))
            ind2=np.argmin(np.array(Qps[100:200]))
            ind3=np.argmin(np.array(Qps[200:300]))
            ind4=np.argmin(np.array(Qps[300:400]))
        else:
            ind1=np.argmax(np.array(Qps[0:100]))
            ind2=np.argmax(np.array(Qps[100:200]))
            ind3=np.argmax(np.array(Qps[200:300]))
            ind4=np.argmax(np.array(Qps[300:400]))
        # quotas
        Qp1=Qps[0:100][ind1]
        Qp2=Qps[100:200][ind2]
        Qp3=Qps[200:300][ind3]
        Qp4=Qps[300:400][ind4]

        # burst size
        bs1=betas[0:100][ind1] 
        bs2=betas[100:200][ind2]
        bs3=betas[200:300][ind3]
        bs4=betas[300:400][ind4]
        betas=[bs1,bs2,bs3,bs4]

        # latent period
        lp1=round(lps[ind1],2)
        lp2=round(lps[100:200][ind2],2)
        lp3=round(lps[200:300][ind3],2)
        lp4=round(lps[300:400][ind4],2)

        lpes=[lp1,lp2,lp3,lp4]

        Qps=[Qp1,Qp2,Qp3,Qp4]

        # maximum growth rates
        mu1=mu_max[0:100][ind1]
        mu2=mu_max[100:200][ind2]
        mu3=mu_max[200:300][ind3]
        mu4=mu_max[300:400][ind4]
        mu_max=[mu1,mu2,mu3,mu4]

        # nutrient half saturation constant
        Nc1=Ncs[0:100][ind1]
        Nc2=Ncs[100:200][ind2]
        Nc3=Ncs[200:300][ind3]
        Nc4=Ncs[300:400][ind4]
        Ncs=[Nc1, Nc2, Nc3, Nc4]

        #volumes
        V1=Vols[0:100][ind1]
        V2=Vols[100:200][ind2]
        V3=Vols[200:300][ind3]
        V4=Vols[300:400][ind4]
        Vs=[V1, V2, V3, V4]
    print(mu_max)

    # choose phytoplankton type
    typePhytos=['Diatom', 'Eukaryote', 'Synechococcus', 'Prochlorochoccus']
    typePhyto=typePhytos[indice]

    bs=betas[indice]

    # temperature dependency
    Temp=20
    T_dep=1 # no T_dep in idealized case 
    suffix=typePhyto+'_BS'+str(round(bs, -1))+'_LOI'+eff
    if m2!=0:
        suffix+='_m2-'+str(m2)
        type_mod='SIVZ'
    else:
        type_mod='SIVZ'
    if dz2==0:
        suffix+='_no-dz2'
    if alph!=1:
        suffix+='_IG'+str(alph)
    if big=='1':
        suffix+='_big'
    if otype != '0':
        suffix+='_'+otype
        tauT = 0.8
        bT=4
        AT=4000
        BT=0.0003
        TN=293.15
        if otype=='upwelling':
            Temp=10
            T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
            R=5
            mu_max=[mu_max[i]*T_dep*R/(R+Ncs[i]) for i in range(N)]
            print(mu_max)
        elif otype=='oligotrophic':
            Temp=25
            T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
            R=0.1
            mu_max=[mu_max[i]*R/(R+Ncs[i])*T_dep for i in range(N)]
            print(mu_max)
        elif otype=='mesotrophic':
            Temp=20
            T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
            R=1
            mu_max=[mu_max[i]*R/(R+Ncs[i])*T_dep for i in range(N)]
            print(mu_max)
    
    # theoretical virus phi
    virus_radius=r_virus[indice]
    Host_volume=Vs[indice]
    phi_e, phi_e_a=phivs_encounter(virus_radius, Host_volume, Temp,indice)

    
    # simulation params
    dt=1/(48)
    nyears=20
    ndays=round(365*nyears)

    mu=mu_max[indice]
    Nc=Ncs[indice]
    mui=0 # infected phyto max growth rate
    beta=float(bs) # burst size
    d=0.1*T_dep # phyto mortality rate
    m=0.1*T_dep #virus mortality
    m2=m2*T_dep #virus qudratic mortality
    Qv=Qvs[indice]
    Qp=Qps[indice]
    eps=1 # adsorption efficiency coefficient
    epso=float(eff) # loss of infectivity rate
    #eps_r=1e-6 # conversion to resistant type
    #eps_lr=1e-6 # loss of resitance
    lat_per=lpes[indice]

    
    N_res=10 # nutrient concentration below MLD
    dN=0.1 # surface-deep mixing rate
    R=0.5 # nutrient concentration
    if otype=='upwelling':
        dN=0.5
        SN=N_res*dN
    elif otype=='mesotrophic':
        dN=0.05
        SN=N_res*dN
    elif otype=='oligotrophic':
        dN=0.01
        SN=N_res*dN
        
    # carrying capacity
    KC_s=(-dN*R+SN)*(R+Nc)/(mu*R)
    CC=KC_s/(mu-d)
    
    print('Carying capacity')
    print(KC_s)
    
    print('Nc')
    print(Nc)
    # grazing parameters
    S_dep=1
    Sc=0.226 #dutkiewicz 2020 (1.5 in umolC.L-1)
    if otype=='upwelling':
        #S_l=5
        S_dep=pow(KC_s,2)/(pow(KC_s,2)+pow(Sc, 2))
    elif otype=='mesotrophic':
        #S_l=1
        S_dep=pow(KC_s,2)/(pow(KC_s,2)+pow(Sc, 2))
    elif otype=='oligotrophic':
        #S_l=0.2
        S_dep=pow(KC_s,2)/(pow(KC_s,2)+pow(Sc, 2))
    
    if pred_type=='encounter_model':
        Vol_host=Vs[indice]
        pi=3.14159265359
        rh=pow(3*Vol_host/(4*pi),1/3)
        phiz=phiz_encounter(rh, rz, Temp)
        phiz=phiz*86400*1000/Qz
        print(phiz)
        suffix+='_encounter_model_zoop'
    else:
        phiz=9.8*T_dep*S_dep
    eps_z=0.3
    dz=0.067*T_dep
    dz2=dz2*T_dep

    # phi and latent period parameter space explored
    if indice in [1,2,3] and big=='0':
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    elif indice in [0] and big=='0':
        phis=list(np.logspace(-11, -7, (12-8)*9+1))
    elif indice==1 and big=='1':
        print('ok')
        phis=list(np.logspace(-9, -4, (13-8)*9+1))
        print(phis)
    elif indice in [0,2,3] and big=='1':
        phis=list(np.logspace(-10, -6, (12-8)*9+1))
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]

    grid_test = [(ph, lp) for ph in phis for lp in lps]

    # explore intraguil predation strength
    alphas=[0, 0.1, 0.5, 1]

    #matrix to store data from theory
    final_SI_th=np.zeros( (len(phis),len(lps)) ) # total phyto: P=S+I
    final_S_th=np.zeros( (len(phis),len(lps)) ) # total S
    final_I_th=np.zeros( (len(phis),len(lps)) ) #total I
    final_V_th=np.zeros( (len(phis),len(lps)) ) # total V
    final_Z_th=np.zeros( (len(phis),len(lps)) ) #total Z
    final_perc_inf_th=np.zeros( (len(phis),len(lps)) ) # % inf
    final_ZV_ratio_th=np.zeros( (len(phis),len(lps)) ) # molar Zv ratio
    final_SI_th[:]=np.nan
    final_S_th[:]=np.nan
    final_I_th[:]=np.nan
    final_V_th[:]=np.nan
    final_Z_th[:]=np.nan
    final_perc_inf_th[:]=np.nan
    final_ZV_ratio_th[:]=np.nan

    n_state=5
    limit_val_mat = np.zeros( (len(phis),len(lps)) ) # matrix for limit of theory validity
    # percentages of virus-induced, zooplankton induced and other (respirtion) induced mortality
    mortality_ratios = np.zeros( (len(phis),len(lps)) ) 
    mortality_ratios_Z = np.zeros( (len(phis),len(lps)) )
    mortality_ratios_O = np.zeros( (len(phis),len(lps)) )
    mortality_ratios_th = np.zeros( (len(phis),len(lps)) )
    mortality_ratios_Z_th = np.zeros( (len(phis),len(lps)) )
    mortality_ratios_O_th = np.zeros( (len(phis),len(lps)) )
    mortality_ratios_th[:]=np.nan
    mortality_ratios_Z_th[:]=np.nan
    mortality_ratios_O_th[:]=np.nan

    # theoretical stability of states
    a_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    o_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    f_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    for l in range(n_state):
        a_stable_states_matrix[l][:]=np.nan
        o_stable_states_matrix[l][:]=np.nan
        f_stable_states_matrix[l][:]=np.nan

    #exluction times
    exclusion_times=np.zeros( (len(phis),len(lps)) )

    # lower and upper limit of range concentration and target concentrations
    A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc=concentration_ranges(indice, otype)

    # matrices for 'realistic' concentrations and other constraints
    valid_concentrations=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_bis=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_ter=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_4=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_5=np.zeros( (len(phis),len(lps)) )

    #matrices for distance to target concentrations
    distance_to_target=np.zeros( (len(phis),len(lps)) )
    distance_to_target_A=np.zeros( (len(phis),len(lps)) )
    distance_to_target_V=np.zeros( (len(phis),len(lps)) )
    distance_to_target_Z=np.zeros( (len(phis),len(lps)) )

    # net primary productivity
    npp_matrix= np.zeros( (len(phis),len(lps)) )

    # matrices to store simulation outputs
    final_S = np.zeros( (len(phis),len(lps)) ) # total S
    final_SI= np.zeros( (len(phis),len(lps)) ) # phyto P=S+I
    final_Surv = np.zeros( (len(phis),len(lps)) ) # matrix of final state/regime reached by the system
    final_Surv_alpha = [np.zeros( (len(phis),len(lps)) ) for i in range(len(alphas))] # for different alphas
    final_ZV_ratio = np.zeros( (len(phis),len(lps)) ) # ZV molar rratio
    final_ZV_ratio_count = np.zeros( (len(phis),len(lps)) ) # ZV count ratio
    final_Z=np.zeros( (len(phis),len(lps)) ) # total Z
    final_V=np.zeros( (len(phis),len(lps)) ) # total V
    final_perc_inf= np.zeros( (len(phis),len(lps)) ) # % Infected
    final_I=np.zeros( (len(phis),len(lps)) ) # total I
    final_pers=np.zeros( (len(phis),len(lps)) ) # period of oscillation (Fourier analysis)
    final_mods=np.zeros( (len(phis),len(lps)) ) # modulo of fourier transform
    final_stab=np.empty( (len(phis),len(lps)) ) # stability of system
    final_stab[:]=np.nan
    final_unstab_lps=[]
    final_unstab_phis=[]
    theor_surv_phi=[]
    theor_surv_lp=[]
    theor_surv_phi_10=[]
    theor_surv_lp_10=[]
    init_conditions=[0.001,0.001, 0.001,0.001]

    # for coexistence analysis
    ntot=6
    matrices_effects_v=[np.zeros( (len(phis),len(lps)) ) for i in range(ntot)]
    matrices_effects_z=[np.zeros( (len(phis),len(lps)) ) for i in range(ntot)]
    for i in range(ntot):
        matrices_effects_v[i][:]=np.nan
        matrices_effects_z[i][:]=np.nan

    cases=np.zeros( (len(phis),len(lps)) )
    cases[:]=np.nan

    # for lopp over the parameter space
    to_plot=[]
    for k in grid_test:
        phi=k[0]
        lp=k[1]

        i=np.where(np.array(phis)==phi)[0][0]
        j=np.where(np.array(lps)==k[1])[0][0]

        # stability theory
        eta=lp
        if dz2!=0:
            if m2==0:
                a_stable, o_stable, f_stable,inits=check_stable_states_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, n_state)
            else:
                a_stable, o_stable, f_stable,inits=check_stable_states_SIVZ_m2(mu, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, n_state)

            for l in range(n_state):
                a_stable_states_matrix[l][i,j]=a_stable[l]
                o_stable_states_matrix[l][i,j]=o_stable[l]
                f_stable_states_matrix[l][i,j]=f_stable[l]

        # simulation
        result=simulation_SIVZ_rk4(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, alph,dt, ndays, init_conditions)
        S=result[0]
        I=result[1]
        V=result[2]
        Z=result[3]
        mZs=result[4]
        mVs=result[5]
        npps=result[8]

        # limit of theory validity
        if phi*KC_s/Qp<10*d and dz2!=0 and m2==0:
            limit_val_mat[i,j]=1

        # theory
        if dz2!=0:
            g=phiz
            a=mu-d+g*(dz/dz2)
            b=1/CC+g*g*eps_z/dz2
            if eff=='0':
                if m2==0:
                    S_star, I_star, V_star, Z_star=equilibrium_SIVZ(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC)
                else:
                    S_star, I_star, V_star, Z_star, case=equilibrium_SIVZ_m2(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,  eps, epso,phiz,eps_z,dz,dz2, CC, 1000)
                    cases[i,j]=case
                    if case==2 and lp<=7:
                        limit_val_mat[i,j]=1
                    elif case==1 and lp<=7:
                        limit_val_mat[i,j]=1
                if S_star>0 and I_star>0 and V_star>0 and Z_star>0:
                    # fill theoretical matrices
                        final_SI_th[i,j]=mt.log10(S_star+I_star)
                        final_S_th[i,j]=mt.log10(S_star)
                        final_I_th[i,j]=mt.log10(I_star)
                        final_V_th[i,j]=mt.log10(V_star)
                        final_Z_th[i,j]=mt.log10(Z_star)
                        final_perc_inf_th[i,j]=I_star*100/(I_star+S_star)
                        final_ZV_ratio_th[i,j]=mt.log10(V_star/Z_star)
                        
                        mZ=phiz*Z_star
                        pi=I_star/(I_star+S_star)
                        mV=pi/eta
                        mortality_ratios_th[i,j]=mV*100/(mZ+mV+d)
                        mortality_ratios_Z_th[i,j]=mZ*100/(mZ+mV+d)
                        mortality_ratios_O_th[i,j]=d*100/(mZ+mV+d)

                if V_star>0 and Z_star>0:
                    theor_surv_phi.append(i)
                    theor_surv_lp.append(len(lps)-j-1)
                if V_star>0 and Z_star>0 and V_star>0:
                    theor_surv_phi_10.append(i)
                    theor_surv_lp_10.append(len(lps)-j-1)
                if V_star>0 and Z_star>0:
                    eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
                    st=0
                    for ei in eigs:
                        ei_r=ei.real
                        if ei_r>0:
                            st=1
                    final_stab[i,j]=st
                    if st==1:
                        final_unstab_phis.append(i)
                        final_unstab_lps.append(len(lps)-j-1)
                if m2==0:
                    for l, alpha in enumerate(alphas):
                        S_star, I_star, V_star, Z_star=equilibrium_SIVZ_alpha(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, alpha)

                        if V_star>0 and Z_star>0:
                            final_Surv_alpha[l][i,j]=1
            else:
                S_star, I_star, V_star, Z_star=equilibrium_SIVZ_LOI(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, eff)
                if V_star>0 and Z_star>0:
                    theor_surv_phi.append(i)
                    theor_surv_lp.append(len(lps)-j-1)
                if V_star>0 and Z_star>0 and V_star>0:
                    theor_surv_phi_10.append(i)
                    theor_surv_lp_10.append(len(lps)-j-1)
                if V_star>0 and Z_star>0:
                    eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
                    st=0
                    for ei in eigs:
                        ei_r=ei.real
                        if ei_r>0:
                            st=1
                    final_stab[i,j]=st
                    if st==1:
                        final_unstab_phis.append(i)
                        final_unstab_lps.append(len(lps)-j-1)
                for l, alpha in enumerate(alphas):
                    S_star, I_star, V_star, Z_star=equilibrium_SIVZ_LOI_alpha(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, eff, alpha)
                    if V_star>Qv and Z_star>Qz:
                        final_Surv_alpha[l][i,j]=1

        
        # extract features of last year of simualtion
        i1=max([0, int(365*(nyears-1)/dt)+1])
        i2=len(V)
        V = list(map(lambda x: 0 if x == float('inf') else x, V))
        I = list(map(lambda x: 0 if x == float('inf') else x, I))
        S = list(map(lambda x: 0 if x == float('inf') else x, S))
        V = list(map(lambda x: 0 if x == float('inf') else x, V))
        I = list(map(lambda x: 0 if x == float('inf') else x, I))
        S = list(map(lambda x: 0 if x == float('inf') else x, S))
        Z = list(map(lambda x: 0 if x == float('inf') else x, Z))
        Sa=np.array(S[i1:i2])
        Ia=np.array(I[i1:i2])
        Va=np.array(V[i1:i2])
        Za=np.array(Z[i1:i2])
        NPPa=np.array(npps[i1:i2])

        # existence criteria
        mZa=np.array(mZs[i1:i2])
        mVa=np.array(mVs[i1:i2])
        surv_S=np.max(Sa/Qp)>1
        surv_V = np.max(Va/Qv)>1 
        surv_Z=np.max(Za/Qz)>1 
        surv_I=np.max(Ia/Qp)>1

        
        # fill matrices
        if surv_S or surv_V or surv_I:
            final_S[i,j]=np.log10(np.nanmean(Sa))
            final_SI[i,j]=np.log10(np.nanmean(Sa+Ia))
            npp_matrix[i,j]=np.log10(np.nanmean(NPPa*106/16))
        else:
            final_S[i,j]=float("NaN")
            final_SI[i,j]=float("NaN")
            npp_matrix[i,j]=float("NaN")
        if surv_V:
            f0, modulos=fft_data(Va,dt)
            mx_m=np.max(modulos)
            if mt.isnan(mx_m):
                mx_per=float("NaN")
            else:
                mx_per=1/f0[np.nanargmax(modulos)]
            final_mods[i,j]=mt.log10(mx_m)
            final_pers[i,j]=mt.log10(mx_per)
        else:
            final_mods[i,j]=float("NaN")
            final_pers[i,j]=float("NaN")

        if not (surv_V or surv_I)  and not surv_Z  and not surv_S:
            surv=0
        elif not (surv_V or surv_I)  and not surv_Z  and surv_S:
            surv=1
        elif (surv_V or surv_I) and not surv_Z:
            surv=2
        elif not (surv_V or surv_I) and surv_Z:
            surv=3
        elif (surv_V or surv_I) and surv_Z:
            surv=4

        # coexistence fitness analysis in the case of coexistence
        if surv==4:
            i1=round(365*(nyears-2)/dt)
            i2=len(V)
            Sf=S[i1:i2]
            If=I[i1:i2]
            Vf=V[i1:i2]
            Zf=Z[i1:i2]
            result_c=coexistence_analysis_SIVZ(Sf, If, Vf, Zf,mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, ntot)

            for l in range(ntot):
                if l==0:
                    matrices_effects_z[l][i,j]=sum(result_c[0:(ntot-1)])
                    matrices_effects_v[l][i,j]=sum(result_c[(ntot-1):(2*ntot-3)])
                if l>0:
                    matrices_effects_z[l][i,j]=result_c[l-1]
                    matrices_effects_v[l][i,j]=result_c[ntot+l-2]
                    
        # exclusion times
        if surv==4 and dz2==0:
            to_plot.append(k)
        t_e=float("NaN")
        if surv==2:
            Zss=np.array(Z)
            Zss=Zss[0::(48*4)]
            for indi in range(len(Zss)):
                if sum(Zss[indi:]/Qz>1)==0:
                    t_e=indi*dt*48*4/365
                    break
        elif surv==3:
            Vss=np.array(V)
            Vss=Vss[0::(48*4)]
            for indi in range(len(Vss)):
                if sum(Vss[indi:]/Qv>1)==0:
                    t_e=indi*dt*48*4/365
                    break
        elif surv in [0,1]:
            Vss=np.array(V)
            Vss=Vss[0::(48*4)]
            Zss=np.array(Z)
            Zss=Zss[0::(48*4)]
            for indi in range(len(Zss)):
                if sum(Zss[indi:]/Qz>1)==0 or sum(Vss[indi:]/Qv>1)==0:
                    t_e=indi*dt*48*4/365
                    break
        exclusion_times[i,j]=t_e

        # fill matrices
        if surv==4:
            mVaZa=np.nanmean(Va/Za)
            final_ZV_ratio[i,j]=np.log10(mVaZa)

        else:
            final_ZV_ratio[i,j]=float("NaN")

        if surv==4:
            mVaZa=np.mean(Va*Qz/(Za*Qv))
            final_ZV_ratio_count[i,j]=np.log10(mVaZa)

            mortality_ratios[i,j]=np.mean(mVa*100/(mZa+mVa+d))
            mortality_ratios_Z[i,j]=np.mean(mZa*100/(mZa+mVa+d))
            mortality_ratios_O[i,j]=np.mean(d*100/(mZa+mVa+d))
        elif surv==2:
            final_ZV_ratio_count[i,j]=float("NaN")
            mortality_ratios[i,j]=np.mean(mVa*100/(mVa+d))
            mortality_ratios_Z[i,j]=float("NaN")
            mortality_ratios_O[i,j]=np.mean(d*100/(mZa+mVa+d))
        elif surv==3:
            final_ZV_ratio_count[i,j]=float("NaN")
            mortality_ratios[i,j]=float("NaN")
            mortality_ratios_Z[i,j]=np.mean(mZa*100/(mZa+mVa+d))
            mortality_ratios_O[i,j]=np.mean(d*100/(mZa+mVa+d))
        else:
            final_ZV_ratio_count[i,j]=float("NaN")
            mortality_ratios[i,j]=float("NaN")
            mortality_ratios_Z[i,j]=float("NaN")
            mortality_ratios_O[i,j]=float("NaN")

        if surv==4 or surv==3:
            final_Z[i,j] = np.log10(np.nanmean(Za))
        else:
            final_Z[i,j]= float("NaN")
        if surv_V or surv_I:
            mVa=np.mean(Va)
            final_V[i,j] = np.log10(mVa)
        else:
            final_V[i,j] = float("NaN")
        if surv_I or surv_V:
            final_I[i,j] = np.log10(np.nanmean(Ia))
        else:
            final_I[i,j] = float("NaN")
        if surv_I or surv_V:
            final_perc_inf[i,j]= np.mean(Ia*100/(Ia+Sa))
        else:
            final_perc_inf[i,j]= float("NaN")
        final_Surv[i,j]=surv

        mean_A=np.nanmean(Sa+Ia)/Qp
        mean_V=np.nanmean(Va)/Qv
        mean_Z=np.nanmean(Za)/Qz
        vals_r=[mean_A, mean_V, mean_Z]

        # check 'realistic' concentrations and distance to target
        if surv==4:
            if mean_A > A_cond_low and mean_A < A_cond_high and mean_V>V_cond_low and mean_V<V_cond_high and mean_Z>Z_cond_low and mean_Z<Z_cond_high:
                valid_concentrations[i,j]=1
            if valid_concentrations[i,j]==1 and final_perc_inf[i,j]>I_cond_low and final_perc_inf[i,j]<I_cond_high:
                valid_concentrations_bis[i,j]=1
            if valid_concentrations[i,j]==1 and mortality_ratios[i,j]>perc_cond_low and mortality_ratios[i,j]<perc_cond_high:
                valid_concentrations_ter[i,j]=1
            if valid_concentrations[i,j]==1 and mean_V/mean_A>=1:
                valid_concentrations_4[i,j]=1
            if valid_concentrations_bis[i,j]==1 and valid_concentrations_ter[i,j]==1 and valid_concentrations_4[i,j]==1:
                valid_concentrations_5[i,j]=1
        
            distance_to_t =  absolute_error(target_conc, vals_r)
            distance_to_target[i,j]=distance_to_t
        else:
            distance_to_t = float("NaN")
            distance_to_target[i,j]=distance_to_t
        if surv_S or surv_V or surv_I:
            distance_to_target_A[i,j]=vals_r[0]-target_conc[0]
        else:
            distance_to_target_A[i,j]=float("NaN")
        if surv_V or surv_I:
            distance_to_target_V[i,j]=vals_r[1]-target_conc[1]
        else:
            distance_to_target_V[i,j]=float("NaN")
        if surv_Z:
            distance_to_target_Z[i,j]=vals_r[2]-target_conc[2]
        else:
            distance_to_target_Z[i,j]=float("NaN")

    # time series of 20 years in the case of the SIVZ model with dv2=dz2=0
    if dz2==0 and m2==0 and nyears==20:
        pp = PdfPages(type_mod+'_model_phi_latent_period_time_series_coex_'+suffix+'.pdf')
        peaks_to_save=[]
        for k in to_plot[0::2]:
            tit='phi='+str(phi)+' lp='+str(lp)
            phi=k[0]
            lp=k[1]
            d=0.1
            nyears0=20
            ndays=365*nyears0
            result=simulation_SIVZ_rk4(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays , init_conditions )
            S=result[0]
            I=result[1]
            V=result[2]
            Z=result[3]
            i1=0
            i2=len(S)-1
            make_1_plot_SIVZ(S, I, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz, 'no')
            i1=0
            i2=round(365*5/dt)
            make_1_plot_SIVZ(S, I, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz, 'no')
            i1=round(365*(nyears0-1)/dt)
            i2=len(S)-1
            make_1_plot_SIVZ(S, I, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz, 'no')
            S_sub=np.array(S[i1:i2])
            I_sub=np.array(I[i1:i2])
            V_sub=np.array(V[i1:i2])
            Z_sub=np.array(Z[i1:i2])

            S_peak, S_peak_l, S_peak_b=find_low_and_high(S_sub)
            I_peak, I_peak_l, I_peak_b=find_low_and_high(I_sub)
            V_peak, V_peak_l, V_peak_b=find_low_and_high(V_sub)
            Z_peak, Z_peak_l, Z_peak_b=find_low_and_high(Z_sub)
            peaks_to_save.append([phi, lp, S_peak*dt, S_peak_l*dt, I_peak*dt, I_peak_l*dt, V_peak*dt, V_peak_l*dt, Z_peak*dt, Z_peak_l*dt])
       
            if not np.isnan(S_peak) and not np.isnan(S_peak_b):
                make_1_plot_SIVZ(S_sub, I_sub, V_sub, Z_sub, int(S_peak), int(S_peak_b), tit,dt,pp, Qp, Qv, Qz, 'yes')
        peaks_to_save=np.array(peaks_to_save)
        write_matrix(peaks_to_save, 'peaks_'+suffix+'.txt', sep=' ')
            
        pp.close()
        
    # time series for specific adsorption rates
    pp = PdfPages(type_mod+'_model_phi_latent_period_time_series_'+suffix+'.pdf')
    if indice in [1,2,3] and big=='0':
        phi_to_plot=[1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    elif indice in [0] and big=='0':
        phi_to_plot=[1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
    elif indice in [0,2,3] and big=='1':
        phi_to_plot=[1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
    elif indice==1 and big=='1':
        phi_to_plot=[1e-9,1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
    for phi in phi_to_plot:
        tit='phi='+str(phi)
        d=0.1
        
        nyears0=5
        ndays=365*nyears0
        result=simulation_SIVZ_rk4(mu, mui, lat_per, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays , init_conditions )
        S=result[0]
        I=result[1]
        V=result[2]
        Z=result[3]
        i1=0
        i2=len(S)-1
        make_plots_SIVZ(S, I, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
        i3=round((365)/dt)
        make_plots_SIVZ(S, I, V, Z, i1, i3,tit,dt,pp, Qp, Qv, Qz)
    pp.close()


    if indice==1 and big=='1':
        lenx=6
    else:
        lenx=5

    atickx=[i*9 if i>0 else 0 for i in range(lenx) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]

    final_S[final_S==-inf]=float("NaN")
    final_Z[final_Z==-inf]=float("NaN")
    final_V[final_V==-inf]=float("NaN")
    final_I[final_I==-inf]=float("NaN")
    final_perc_inf[final_perc_inf==-inf]=float("NaN")

    # pdf with figures of matrices exploring the parameter space
    pp = PdfPages(type_mod+'_model_phi_latent_period_'+suffix+'.pdf')
    ylb='Latent period'
    bounds=[0,1,2,3,4]
    norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)
    mi=0
    mx=4
    n_state_R=11
    mx_R=n_state_R-1
    colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx_R+2)])
    cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
    colors.append(cn)
    idis=[0,2,6,4,9]
    colors=[colors[idi] for idi in idis]
    colors=tuple(colors)
    cmap_R= matplotlib.colors.ListedColormap(colors)
    axi=plot_with_scale(final_Surv,cmap_R,mi,mx, atickx, aticky, alabel_tickx, alabel_ticky, 'State reached', norm=norm, yl=ylb)
    coex_mat=np.zeros( (len(phis),len(lps)) )
    coex_mat[final_Surv==4]=1
    coex_mat=np.transpose(coex_mat)
    coex_mat=np.flipud(coex_mat)
    draw_borders(coex_mat, 1.5, 0.05, 'white', axi)
    pp.savefig()

    # theory (when dz2!=0, dv2==0 or dv2!=0)
    if dz2!=0:
        mi=0
        mx=2
        cmap = matplotlib.colors.ListedColormap(['red', 'green', 'dodgerblue'])
        eq_types=['S,I,V,Z', 'S,I,V,0', 'S,0,0,Z', 'S,0,0,0', '0,0,0,0']
        for l in range(n_state):
            matri=copy.deepcopy(f_stable_states_matrix[l])
            matri[limit_val_mat==0]=np.nan
            plot_with_scale_bis(matri, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, eq_types[l], yl=ylb)
            pp.savefig()
        predicted_state=np.empty((len(phis), len(lps)))
        predicted_state[:]=np.nan
        jo=n_state-1
        to_check=list(range(n_state-1))
        for l in range(n_state):
            tab=f_stable_states_matrix[l]
            c_lower=[]
            c_lower_bis=[]
            if l in to_check: # check for alternate stable state
                for j in range(1, n_state-l):
                    tab_min1=f_stable_states_matrix[l+j]
                    condi = tab_min1!=2
                    condi=condi*1
                    c_lower.append(condi)

                    condi = tab_min1!=1
                    condi=condi*1
                    c_lower_bis.append(condi)
            c1=tab==2
            c1_i=c1*1
            c2=np.isnan(predicted_state)
            c2_i=c2*1
            c3=np.array(c2_i*c1_i, dtype=bool)
            predicted_state[c3] = jo

            c1=tab==1
            c1_i=c1*1
            c2=np.isnan(predicted_state)
            c2_i=c2*1
            c3=np.array(c2_i*c1_i, dtype=bool)
            predicted_state[c3] = jo

            c1=tab==0
            c1_i=c1*1
            c2=np.isnan(predicted_state)
            c2_i=c2*1
            if l in to_check:
                for j in range(len(c_lower)):
                    c2_i=c2_i*c_lower[j]*c_lower_bis[j]
            c3=np.array(c2_i*c1_i, dtype=bool)
            predicted_state[c3] = jo
            jo=jo-1

        alter_state=np.zeros((len(phis), len(lps)))
        u_state=np.zeros((len(phis), len(lps)))
        j0=n_state-1
        for k in range(n_state):
            tabi=f_stable_states_matrix[k]
            tabi[np.isnan(tabi)]=0
            condi = tabi==2
            condi_bis= tabi==1
            condi_ter=np.logical_or(condi, condi_bis)
            condi_ter=condi_ter*1
            alter_state+=condi_ter

            condi = tabi==0
            condi=condi*1
            u_state+=condi
            j0=j0-1

        alter_state[limit_val_mat==0]=0
        u_state[limit_val_mat==0]=0

        alter_state[alter_state==1]=0
        alter_state[alter_state>=1]=1

        u_state[u_state!=n_state]=0
        u_state[u_state==n_state]=1

        # convert to the state from simulation
        predicted_state_bis=np.zeros((len(phis), len(lps)))
        predicted_state_bis[:]=np.nan
        for j in range(n_state):
            if j==0:
                predicted_state_bis[predicted_state==j]=0
            if j==1:
                predicted_state_bis[predicted_state==j]=1    
            if j==2:
                predicted_state_bis[predicted_state==j]=3
            if j==3:
                predicted_state_bis[predicted_state==j]=2
            if j==4:
                predicted_state_bis[predicted_state==j]=4

        mi=0
        mx=4       
        alter_state_bis=copy.deepcopy(alter_state)
        alter_state_bis[alter_state_bis==1]=2
        alter_state_bis[alter_state_bis==0]=1
        alter_state_bis[alter_state_bis==2]=0
        predicted_state_bis[limit_val_mat==0]=np.nan
        axi=plot_with_scale(predicted_state_bis,cmap_R ,mi,mx,atickx, aticky, alabel_tickx, alabel_ticky, 'Predicted state', norm=norm, yl=ylb)
        plt.contour(np.transpose(alter_state), colors=['limegreen', 'white'], levels=[0,1] , extent=[-1, len(phis), -1, len(lps)], origin='upper', linewidths=0.5)
        u_state_bis=copy.deepcopy(u_state)
        u_state_bis[u_state_bis==1]=2
        u_state_bis[u_state_bis==0]=1
        u_state_bis[u_state_bis==2]=0
        plt.contourf(np.transpose(alter_state_bis), colors=['limegreen', 'white'], levels=[0,0.9] , alpha=0.2, extent=[-1, len(phis), -1, len(lps)], origin='upper')
        plt.contour(np.transpose(u_state), colors=['black', 'white'], levels=[0,1] , extent=[-1, len(phis), -1, len(lps)], origin='upper', linewidths=0.5)
        plt.contourf(np.transpose(u_state_bis), colors=['black', 'white'], levels=[0,0.9] , alpha=0.2, extent=[-1, len(phis), -1, len(lps)], origin='upper')
        
        coex_mat=np.zeros( (len(phis),len(lps)) )
        coex_mat[predicted_state_bis==4]=1
        coex_mat=np.transpose(coex_mat)
        coex_mat=np.flipud(coex_mat)
        draw_borders(coex_mat, 1.5, 0.05, 'white', axi) 
        pp.savefig()

    mi=np.nanmin(np.array(exclusion_times))
    mx=np.nanmax(np.array(exclusion_times))
    plot_with_scale(exclusion_times, 'YlGnBu', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Exclusion time (year)', yl=ylb)
    pp.savefig()

    exclusion_times_bis=copy.deepcopy(exclusion_times)
    exclusion_times_bis[exclusion_times_bis>1.5]=1.5
    plot_with_scale(exclusion_times_bis, 'YlGnBu', mi, 1.5, atickx, aticky, alabel_tickx, alabel_ticky , 'Exclusion time (year)', yl=ylb)
    pp.savefig()

    mi=0
    mx=100
    plot_with_scale_bis(mortality_ratios, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% virus induced death', yl=ylb)
    pp.savefig()
    plot_with_scale_bis(mortality_ratios_Z, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% zoop induced death', yl=ylb)
    pp.savefig()
    plot_with_scale_bis(mortality_ratios_O, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% other induced death', yl=ylb)
    pp.savefig()

    mortality_ratios_th[limit_val_mat==0]=np.nan
    mortality_ratios_Z_th[limit_val_mat==0]=np.nan
    mortality_ratios_O_th[limit_val_mat==0]=np.nan
    plot_with_scale_bis(mortality_ratios_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% virus induced death th', yl=ylb)
    pp.savefig()
    plot_with_scale_bis(mortality_ratios_Z_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% zoop induced death th', yl=ylb)
    pp.savefig()
    plot_with_scale_bis(mortality_ratios_O_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% other induced death th', yl=ylb)
    pp.savefig()


    absmax=max([abs(np.nanmin(np.array(final_ZV_ratio))), abs(np.nanmax(np.array(final_ZV_ratio)))])
    mi=-absmax
    mx=absmax
    cmap=matplotlib.colormaps['bwr']
    rcmap = cmap.reversed()
    plot_with_scale_bis(final_ZV_ratio, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Zooplankton)', yl=ylb)
    pp.savefig()

    absmax=max([abs(np.nanmin(np.array(final_ZV_ratio_count))), abs(np.nanmax(np.array(final_ZV_ratio_count)))])
    mi=-absmax
    mx=absmax
    plot_with_scale_bis(final_ZV_ratio_count, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Zooplankton) count', yl=ylb)
    pp.savefig()

    final_ZV_ratio_th[limit_val_mat==0]=np.nan
    absmax=max([abs(np.nanmin(np.array(final_ZV_ratio_th))), abs(np.nanmax(np.array(final_ZV_ratio_th)))])
    mi=-absmax
    mx=absmax
    plot_with_scale_bis(final_ZV_ratio_th, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Zooplankton) th', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z))
    mx=np.nanmax(np.array(final_Z))
    plot_with_scale(final_Z, 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z)-mt.log10(Qz))
    mx=np.nanmax(np.array(final_Z)-mt.log10(Qz))
    plot_with_scale(np.array(final_Z)-mt.log10(Qz), 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    final_Z_th[limit_val_mat==0]=np.nan
    mi=np.nanmin(np.array(final_Z_th)-mt.log10(Qz))
    mx=np.nanmax(np.array(final_Z_th)-mt.log10(Qz))
    plot_with_scale(np.array(final_Z_th)-mt.log10(Qz), 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(ind.L-1)) th', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V))
    mx=np.nanmax(np.array(final_V))
    plot_with_scale(final_V, 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V)-mt.log10(Qv))
    mx=np.nanmax(np.array(final_V)-mt.log10(Qv))
    plot_with_scale(np.array(final_V)-mt.log10(Qv), 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    final_V_th[limit_val_mat==0]=np.nan
    mi=np.nanmin(np.array(final_V_th)-mt.log10(Qv))
    mx=np.nanmax(np.array(final_V_th)-mt.log10(Qv))
    plot_with_scale(np.array(final_V_th)-mt.log10(Qv), 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus th (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_SI))
    mx=np.nanmax(np.array(final_SI))
    plot_with_scale(final_SI, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_SI)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_SI)-mt.log10(Qp))
    plot_with_scale(np.array(final_SI)-mt.log10(Qp), 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    final_SI_th[limit_val_mat==0]=np.nan
    mi=np.nanmin(np.array(final_SI_th)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_SI_th)-mt.log10(Qp))
    plot_with_scale(np.array(final_SI_th)-mt.log10(Qp), 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I th (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_S))
    mx=np.nanmax(np.array(final_S))
    plot_with_scale(final_S, 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_S)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_S)-mt.log10(Qp))
    plot_with_scale(np.array(final_S)-mt.log10(Qp), 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    final_S_th[limit_val_mat==0]=np.nan
    mi=np.nanmin(np.array(final_S_th)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_S_th)-mt.log10(Qp))
    plot_with_scale(np.array(final_S_th)-mt.log10(Qp), 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible th (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_I))
    mx=np.nanmax(np.array(final_I))
    plot_with_scale(final_I, 'Oranges', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Infected (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_I)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_I)-mt.log10(Qp))
    plot_with_scale(np.array(final_I)-mt.log10(Qp), 'Oranges', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Infected (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    final_I_th[limit_val_mat==0]=np.nan
    mi=np.nanmin(np.array(final_I_th)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_I_th)-mt.log10(Qp))
    plot_with_scale(np.array(final_I_th)-mt.log10(Qp), 'Oranges', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Infected th (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=0
    mx=100
    plot_with_scale_bis(final_perc_inf, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , '% of infected cells', yl=ylb)
    pp.savefig()

    final_perc_inf_th[limit_val_mat==0]=np.nan
    plot_with_scale_bis(final_perc_inf_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , '% of infected cells th', yl=ylb)
    pp.savefig()

    
    mi=np.nanmin(np.array(npp_matrix))
    mx=np.nanmax(np.array(npp_matrix))
    plot_with_scale(np.array(npp_matrix), 'YlGn', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'NPP (log10(umolC.L-1.d-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_pers))
    mx=np.nanmax(np.array(final_pers))
    plot_with_scale(final_pers, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Dominant period (log10(days))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_mods))
    mx=np.nanmax(np.array(final_mods))
    plot_with_scale(final_mods, 'inferno', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'FFT modulo', yl=ylb)
    pp.savefig()

    mi=0
    mx=1
    plot_with_scale(cases, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Approx type', yl=ylb)
    pp.savefig()

    mi=0
    mx=400
    distance_to_target_bis=copy.deepcopy(distance_to_target)
    distance_to_target_bis[distance_to_target_bis>mx]=mx
    plot_with_scale_bis(distance_to_target_bis, 'inferno', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target concentrations', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(distance_to_target_A))
    mx=np.nanmax(np.array(distance_to_target_A))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_A, 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target A', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(distance_to_target_V))
    mx=np.nanmax(np.array(distance_to_target_V))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_V, 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target V', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(distance_to_target_Z))
    mx=np.nanmax(np.array(distance_to_target_Z))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_Z, 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target Z', yl=ylb)
    pp.savefig()
    

    mi=np.nanmin(np.array(distance_to_target_A*100/target_conc[0]))
    mx=np.nanmax(np.array(distance_to_target_A*100/target_conc[0]))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_A*100/target_conc[0], 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target A (%)', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(distance_to_target_V*100/target_conc[1]))
    mx=np.nanmax(np.array(distance_to_target_V*100/target_conc[1]))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_V*100/target_conc[1], 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target V (%)', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(distance_to_target_Z*100/target_conc[2]))
    mx=np.nanmax(np.array(distance_to_target_Z*100/target_conc[2]))
    amx=max([abs(mi),abs(mx)])
    plot_with_scale_bis(distance_to_target_Z*100/target_conc[2], 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, 'Distance to target Z (%)', yl=ylb)
    pp.savefig()


    to_plot=[valid_concentrations, valid_concentrations_bis, valid_concentrations_ter, valid_concentrations_4, valid_concentrations_5]
    titls=['Valid concentrations', 'Valid concentrations+ % infected', 'Valid concentrations + % virus kill','Valid concentrations + V/A ratio' ,'Valid concentrations + % infected + % virus kill + V/A ratio']
    mi=0
    mx=1
    cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue'])
    for k,mato in enumerate(to_plot):
        plot_with_scale(mato, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, titls[k], yl=ylb)
        pp.savefig()

    try:
        ind_min=np.nanargmin(distance_to_target)
    except:
        ind_min, min_d, phi_min,rat_min=np.nan, np.nan, np.nan, np.nan

    if not np.isnan(ind_min):
        ind_min=np.nanargmin(distance_to_target)
        phi_index, lp_index = np.unravel_index(ind_min, distance_to_target.shape)
        min_d=distance_to_target[phi_index, lp_index]
        phi_min=phis[phi_index]
        lp_min=lps[lp_index]
    
        distance_to_target_val=copy.deepcopy(distance_to_target)
        distance_to_target_val[valid_concentrations_5==0]=np.nan
    
        try:
            ind_min_val=np.nanargmin(distance_to_target_val)
        except:
            ind_min_val,min_d_val, phi_min_val, lp_min_val=np.nan,np.nan,np.nan, np.nan
    
        if not np.isnan(ind_min_val):
            phi_index_val, lp_index_val = np.unravel_index(ind_min_val, distance_to_target.shape)
            min_d_val=distance_to_target[phi_index_val, lp_index_val]
            phi_min_val=phis[phi_index_val]
            lp_min_val=lps[lp_index_val]

            valid_concentrations_5[phi_index_val, lp_index_val]=2
            mi=0
            mx=2
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold'])
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()
    
            phi_index_par, lp_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(lat_per-np.array(lps)))
            valid_concentrations_5[phi_index_par, lp_index_par]=3
            mi=0
            mx=3
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 10) & (np.nan_to_num(valid_concentrations_5, nan=np.inf) != 3)
            valid_concentrations_5[cond]=2
            vl=distance_to_target_val[valid_concentrations_5==3]<min_d_val+10
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            
        
            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_5[phi_a_index_par, lp_index_par]=4
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            mx=4
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            if phi_a_index_par==phi_index_par:
                valid_concentrations_5[phi_a_index_par, lp_index_par]=3
            cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) > min_d_val + 10) & (np.nan_to_num(valid_concentrations_5, nan=np.inf) != 3) & (np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 50)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan', 'darkorange'])
            valid_concentrations_5[cond]=5
            mx=5
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

        else:
            phi_index_par, lp_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(lat_per-np.array(lps)))
            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_5[phi_index_par, lp_index_par]=3
            valid_concentrations_5[phi_a_index_par, lp_index_par]=4
            mi=0
            mx=4
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()
            
    else:
        phi_index_par, lp_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(lat_per-np.array(lps)))
        phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
        valid_concentrations_5[phi_index_par, lp_index_par]=3
        valid_concentrations_5[phi_a_index_par, lp_index_par]=4
        mi=0
        mx=4
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
        pp.savefig()

    pp.close()

    if 'vl' not in globals():
        vl=np.nan
    else:
        vl=vl[0]

    to_write=[min_d, phi_min,lp_min, min_d_val, phi_min_val, lp_min_val, phi_e, phi_e_a, lat_per, vl]
    if nyears==20:
        write_vector(to_write, 'optimums_SIVZ_'+suffix+'.txt', ' ')

    # coexistece analysis figures
    pp = PdfPages(type_mod+'_model_phi_latent_period_'+suffix+'_coexistence.pdf')
    atickx=[i*9 if i>0 else 0 for i in range(lenx) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]


    mxs=[]
    mis=[]
    tr=1000
    for l in range(ntot):
        mat1=matrices_effects_z[l]
        mat2=matrices_effects_v[l]
        mat1[mat1<-tr]=-tr
        mat2[mat2<-tr]=-tr
        mat1[mat1>tr]=tr
        mat2[mat2>tr]=tr
        mxs.append(np.max(mat1))
        mxs.append(np.max(mat2))

        mis.append(np.min(mat1))
        mis.append(np.min(mat2))

    mx=max(mxs)
    mi=min(mis)
    amx=max([abs(mi), abs(mx)])
    ma=amx
    nma=-amx

    cmap = matplotlib.colormaps['PuOr_r']
    colors0 = cmap(np.linspace(0, 1, 100))
    cmap0=matplotlib.colors.ListedColormap(colors0)
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)

    titles=['Average growth rate','Fluctuation free growth rate', 'Relative non linearity in S', 'Relative non linearity in I', 'I, S covariance', 'I, S variance interaction']
    titles_save=['Average_growth_rate','Fluctuation_free_growth_rate', 'Relative_linearity_S', 'Relative_non_linearity_I', 'IS_covariance', 'IS_variance_interaction']
    for l in range(ntot):
        mat=matrices_effects_z[l]
        write_matrix(mat, 'coexistence_analysis_Z_'+type_mod+'_'+suffix+'_'+titles_save[l]+'.txt', ' ')
        mat[mat<-tr]=-tr
        mat[mat>tr]=tr
        maxi=np.nanmax(mat)
        mini=np.nanmin(mat)
        amx=max([abs(mini), abs(maxi)])
        ma=amx
        nma=-amx
        if (maxi==mini):
            ma=1
            nma=-1
        bounds=np.linspace(-amx, amx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        plot_with_scale(mat,cmap0,nma,ma, atickx, aticky, alabel_tickx, alabel_ticky, titles[l]+' (Z)', norm=norm, yl=ylb)
        pp.savefig()
    for l in range(ntot):
        mat=matrices_effects_v[l]
        write_matrix(mat, 'coexistence_analysis_V_'+type_mod+'_'+suffix+'_'+titles_save[l]+'.txt', ' ')
        mat[mat<-tr]=-tr
        mat[mat>tr]=tr
        maxi=np.nanmax(mat)
        mini=np.nanmin(mat)
        amx=max([abs(mini), abs(maxi)])
        ma=amx
        nma=-amx
        if (maxi==mini):
            ma=1
            nma=-1
        bounds=np.linspace(-amx, amx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        plot_with_scale(mat,cmap0,nma,ma, atickx, aticky, alabel_tickx, alabel_ticky, titles[l]+' (V)', norm=norm, yl=ylb)
        pp.savefig()
    pp.close()

    # write data to create scaled figures
    if nyears==20:
        write_matrix(exclusion_times, 'exclusion_times_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_S-mt.log10(Qp), 'model_data/Susceptible_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_I-mt.log10(Qp), 'model_data/Infected_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_SI-mt.log10(Qp), 'model_data/Tot_algae_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_V-mt.log10(Qv), 'model_data/Virus_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_Z-mt.log10(Qz), 'model_data/Zoop_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_ZV_ratio, 'model_data/ZV_ratio_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(distance_to_target, 'model_data/distances_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(distance_to_target_A, 'model_data/distance_A_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(distance_to_target_V, 'model_data/distance_V_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(distance_to_target_Z, 'model_data/distance_Z_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(valid_concentrations_5, 'model_data/validation_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_pers, 'model_data/Period_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(final_mods, 'model_data/FFT_modulo_'+type_mod+'_'+suffix+'.txt', ' ')
        write_matrix(npp_matrix, 'model_data/NPP_'+type_mod+'_'+suffix+'.txt', ' ')

    write_matrix(final_pers, 'model_data/Period_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_mods, 'model_data/FFT_modulo_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_S_th-mt.log10(Qp), 'model_data/Susceptible_th_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_I_th-mt.log10(Qp), 'model_data/Infected_th_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_SI_th-mt.log10(Qp), 'model_data/Tot_algae_th_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_V_th-mt.log10(Qv), 'model_data/Virus_th_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_Z_th-mt.log10(Qz), 'model_data/Zoop_th_'+type_mod+'_'+suffix+'.txt', ' ')
    write_matrix(final_ZV_ratio_th, 'model_data/ZV_ratio_th_'+type_mod+'_'+suffix+'.txt', ' ')
    print("=== Script finished ===")
