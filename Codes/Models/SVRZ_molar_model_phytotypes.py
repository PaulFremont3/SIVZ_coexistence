import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import copy
from numpy import inf
from SVRZ_functions import *
from generic_functions import *


if __name__ == '__main__':
    indice=int(sys.argv[1]) # choose phytoplankton type
    otype=sys.argv[2] # ocean type
    gr_ratio=float(sys.argv[3]) # cost of resistance (between 0 and 1)
    dz2=float(sys.argv[4]) # quadratic mortality of zooplankton
    m2=float(sys.argv[5]) # quadratic mortality of virus
    param=str(sys.argv[6]) #resistance type: phir or epsr
    spec_value=str(sys.argv[7]) #specify a specific value of phi/phir or eps/epsr => plot 1d matrix. If 'no' => explore the phi/phir or eps/epsr space
    
    if spec_value!='no':
        spec_value=float(spec_value)

    eff=0
    type_m='SVRZ'
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
    # grazer quota
    Qz=Q_grazer(rz)

    # virus quota
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
    # extract life history traist of smalleest of each type of pytoplankton
    if N==4:
        ind1=np.argmin(np.array(Qps[0:100]))
        ind2=np.argmin(np.array(Qps[100:200]))
        ind3=np.argmin(np.array(Qps[200:300]))
        ind4=np.argmin(np.array(Qps[300:400]))
        Qp1=Qps[0:100][ind1]
        Qp2=Qps[100:200][ind2]
        Qp3=Qps[200:300][ind3]
        Qp4=Qps[300:400][ind4]

        bs1=round(betas[ind1],-1) #
        bs2=round(betas[100:200][ind2],-1)
        bs3=round(betas[200:300][ind3],-1)
        bs4=round(betas[300:400][ind4],-1)
        betas=[bs1,bs2,bs3,bs4]

        lp1=round(lps[ind1],2)
        lp2=round(lps[100:200][ind2],2)
        lp3=round(lps[200:300][ind3],2)
        lp4=round(lps[300:400][ind4],2)

        lpes=[lp1,lp2,lp3,lp4]

        Qps=[Qp1,Qp2,Qp3,Qp4]

        mu1=mu_max[0:100][ind1]
        mu2=mu_max[100:200][ind2]
        mu3=mu_max[200:300][ind3]
        mu4=mu_max[300:400][ind4]
        mu_max=[mu1,mu2,mu3,mu4]

        Nc1=Ncs[0:100][ind1]
        Nc2=Ncs[100:200][ind2]
        Nc3=Ncs[200:300][ind3]
        Nc4=Ncs[300:400][ind4]
        Ncs=[Nc1, Nc2, Nc3, Nc4]

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

    # default temperature
    Temp=20
    T_dep=1 # no T_dep in idealized case
    
    # suffix to identify output pdfs
    suffix=typePhyto+'_BS'+str(bs)+'_LOI'+str(eff)+'_GR-R'+str(gr_ratio)
    if dz2==0:
        suffix+='_no-dz2'
    if m2!=0:
        suffix+='_m2-'+str(m2)

    R=0.5 # nutrient concentration
    # temperature dependencey
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
    mur=mu*float(gr_ratio)
    mui=0 # infected phyto max growth rate
    beta=float(bs) # burst size
    d=0.1*T_dep # phyto mortality rate
    m=0.1*T_dep #virus mortality
    m2=m2*T_dep #virus qudratic mortality
    Qv=Qvs[indice]
    Qp=Qps[indice]
    eps=1 # adsorption efficiency coefficient
    epso=float(eff) # loss of infectivity rate
    eps_r=1e-6 # conversion to resistant type
    eps_lr=1e-6 # loss of resitance
    
    lat_per=lpes[indice]

    N_res=10 # nutrient concentration below MLD
    dN=0.1 # surface-deep mixing rate
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
    KC_s=(-dN*R+SN)/mu
    CC=KC_s/(mu-d)

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
    phiz=9.8*T_dep*S_dep
    eps_z=0.3
    dz=0.067*T_dep
    dz2=dz2*T_dep

    print('Carying capacity')
    print(KC_s)

    # phi parameter sapce and resistance strength parameter space
    if indice in [1,2,3]:
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    elif indice in [0]:
        phis=list(np.logspace(-11, -7, (12-8)*9+1))

    if spec_value=='no':
        phi_over_phir=list(np.logspace(0, 4, 4*9+1))
        phi_over_phir=phi_over_phir[0:len(phi_over_phir)-1]
        phi_over_phir.append(0)
        del phi_over_phir[0]
        print(phi_over_phir)
        suffix+='_'+param
    else:
        phi_over_phir=[spec_value]
        if spec_value==0:
            suffix+='_'+param+'-fullres'
        else:
            suffix+='_'+param+'-'+str(spec_value)

    # parameter space
    grid_test = [(ph, ratio) for ph in phis for ratio in phi_over_phir]

    # matrices to store results (see SIVZ file for specs)
    limit_val_mat = np.zeros( (len(phis),len(phi_over_phir)) )

    mortality_ratios=np.zeros( (len(phis),len(phi_over_phir)) )

    final_state_all=np.empty( (len(phis),len(phi_over_phir)) )
    final_state_all[:]=np.nan

    exclusion_times=np.zeros( (len(phis),len(phi_over_phir)) )

    n_state=11
    a_stable_states_matrix = [np.empty( (len(phis),len(phi_over_phir)) ) for i in range(n_state)]
    o_stable_states_matrix = [np.empty( (len(phis),len(phi_over_phir)) ) for i in range(n_state)]
    f_stable_states_matrix = [np.empty( (len(phis),len(phi_over_phir)) ) for i in range(n_state)]
    for l in range(n_state):
        a_stable_states_matrix[l][:]=np.nan
        o_stable_states_matrix[l][:]=np.nan
        f_stable_states_matrix[l][:]=np.nan


    alphas=[0 , 0.1, 0.5, 1]
    final_S = np.zeros( (len(phis),len(phi_over_phir)) )
    final_Surv = np.zeros( (len(phis),len(phi_over_phir)) )
    final_Surv_alpha = [np.zeros( (len(phis),len(phi_over_phir)) ) for i in range(len(alphas))]
    final_Surv_alpha0 = [np.zeros( (len(phis),len(phi_over_phir)) ) for i in range(len(alphas))]
    final_Surv_alpha1 = [np.zeros( (len(phis),len(phi_over_phir)) ) for i in range(len(alphas))]
    final_ZV_ratio = np.zeros( (len(phis),len(phi_over_phir)) )
    final_ZV_ratio_count = np.zeros( (len(phis),len(phi_over_phir)) )
    final_Z=np.zeros( (len(phis),len(phi_over_phir)) )
    final_V=np.zeros( (len(phis),len(phi_over_phir)) )
    final_I=np.zeros( (len(phis),len(phi_over_phir)) )
    final_perc_inf= np.zeros( (len(phis),len(phi_over_phir)) )
    final_perc_res= np.zeros( (len(phis),len(phi_over_phir)) ) # percentage of resistant cells
    final_R=np.zeros( (len(phis),len(phi_over_phir)) ) # total R
    final_pers=np.zeros( (len(phis),len(phi_over_phir)) )
    final_mods=np.zeros( (len(phis),len(phi_over_phir)) )
    final_stab=np.empty( (len(phis),len(phi_over_phir)) )
    final_stab[:]=np.nan
    final_unstab_phi_over_phirs=[]
    final_unstab_phis=[]
    theor_surv_phi_over_phir=[]
    theor_surv_phi=[]
    theor_surv_phi_all=[]
    theor_surv_phi_over_phir_all=[]
    theor_surv_phi0=[]
    theor_surv_phi_over_phir0=[]
    theor_surv_phi1=[]
    theor_surv_phi_over_phir1=[]

    to_plot=[]

    # realistic and target concentrations
    A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc=concentration_ranges(indice, otype)

    valid_concentrations=np.zeros( (len(phis),len(phi_over_phir)) )
    valid_concentrations_bis=np.zeros( (len(phis),len(phi_over_phir)) )
    valid_concentrations_ter=np.zeros( (len(phis),len(phi_over_phir)) )
    valid_concentrations_4=np.zeros( (len(phis),len(phi_over_phir)) )

    distance_to_target=np.zeros( (len(phis),len(phi_over_phir)) )
    distance_to_target_A=np.zeros( (len(phis),len(phi_over_phir)) )
    distance_to_target_V=np.zeros( (len(phis),len(phi_over_phir)) )
    distance_to_target_Z=np.zeros( (len(phis),len(phi_over_phir)) )

    # main loop through the parameter space
    init_conditions=[0.001,0.001, 0.001, 0.001]
    for k in grid_test:
        phi=k[0]
        if param=='phir':
          if k[1] !=0:
            phir=phi/k[1]
            epsr=eps
          else:
            phir=0
            epsr=eps
        elif param=='epsr':
          if k[1] !=0:
            phir=phi
            epsr=eps/k[1]
          else:
            phir=phi
            epsr=0

        i0=np.where(np.array(phis)==phi)[0][0]
        j0=np.where(np.array(phi_over_phir)==k[1])[0][0]
        alpha=1
        eta=0
        # simulation
        result=simulation_SVRZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur, dt, ndays, init_conditions)
    
        S=result[0]
        V=result[1]
        R=result[2]
        Z=result[3]
        mZ=result[4]
        mV=result[5]

        V = list(map(lambda x: 0 if x == float('inf') else x, V))
        S = list(map(lambda x: 0 if x == float('inf') else x, S))
        R = list(map(lambda x: 0 if x == float('inf') else x, R))
        Z = list(map(lambda x: 0 if x == float('inf') else x, Z))

        # extract last year of simulation
        i1=int(365*(nyears-1)/dt)+1
        i2=len(V)
        Sa=np.array(S[i1:i2])
        Ra=np.array(R[i1:i2])
        Va=np.array(V[i1:i2])
        Za=np.array(Z[i1:i2])

        mZa=np.array(mZ[i1:i2])
        mVa=np.array(mV[i1:i2])    

        # existence criterion
        surv_S=np.max(Sa/Qp)>1 
        surv_R=np.max(Ra/Qp)>1
        surv_V = np.max(Va/Qv)>1 
        surv_Z=np.max(Za/Qz)>1 

        # fill matrices
        if surv_V  or surv_S or surv_R:
            final_S[i0,j0]=np.log10(np.mean(Sa))
            final_perc_res[i0,j0]= np.mean(Ra*100/(Ra+Sa))
        else:
            final_S[i0,j0]=float("NaN")
            final_perc_res[i0,j0]=float("NaN")

        # criteerion to determine if S or R is in sufficient abundance
        surv_S_bis = np.mean(Sa) > 0.00001 
        surv_R_bis = np.mean(Ra) > 0.00001 

        # fourier analysis
        if surv_V:
            f0, modulos=fft_data(Va,dt)
            mx_m=np.max(modulos)
            mx_per=1/f0[np.nanargmax(modulos)]
            final_mods[i0,j0]=mt.log10(mx_m)
            final_pers[i0,j0]=mt.log10(mx_per)
        else:
            final_mods[i0,j0]=float("NaN")
            final_pers[i0,j0]=float("NaN")

        # final state of the system
        surv=np.nan
        if not surv_V  and not surv_Z  and not surv_S and not surv_R:
            final_state_all[i0,j0]=0
            surv=0
        elif surv_R and not surv_V  and not surv_Z  and not surv_S_bis:
            final_state_all[i0,j0]=1
            surv=1
        elif surv_S and not surv_R_bis and not surv_V  and not surv_Z :
            final_state_all[i0,j0]=2
            surv=1
        elif surv_R and not surv_V  and surv_Z and not surv_S_bis:
            final_state_all[i0,j0]=3
            surv=2
        elif surv_S and not surv_R_bis and not surv_V  and surv_Z:
            final_state_all[i0,j0]=4
            surv=2
        elif surv_R and  surv_V  and not surv_Z and not surv_S_bis:
            final_state_all[i0,j0]=5
            surv=3
        elif surv_S and  surv_V  and not surv_Z and not surv_R_bis:
            final_state_all[i0,j0]=6
            surv=3
        elif surv_S and  surv_V  and not surv_Z and surv_R:
            final_state_all[i0,j0]=7
            surv=4
        elif not surv_S_bis and  surv_V  and surv_Z and surv_R:
            final_state_all[i0,j0]=8
            surv=5
        elif surv_S and surv_V  and surv_Z and not surv_R_bis:
            final_state_all[i0,j0]=9
            surv=6
        elif surv_S and surv_V and surv_Z and surv_R:
            final_state_all[i0,j0]=10
            surv=7
        # check if we covered all cases
        if np.isnan(surv):
            surv=0
            final_state_all[i0,j0]=0
            print(surv_S)
            print(surv_V)
            print(surv_Z)
            print(surv_R)
            print(surv_R_bis)
            print(surv_S_bis)

        if surv==7 and dz2==0:
            to_plot.append(k)

        # exclusion time
        t_e=float("NaN")
        if surv in [3,4]:
            Zss=np.array(Z)
            Zss=Zss[0::(48*4)]
            for indi in range(len(Zss)):
                if sum(Zss[indi:]/Qz>1)==0:
                    t_e=indi*dt*48*4/365
                    break
        elif surv ==2:
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
        exclusion_times[i0,j0]=t_e

        # fill matrices
        states_with_VandZ=[7, 6, 5]
        if surv in states_with_VandZ:
            mVaZa=np.nanmean(Va/Za)
            rat=np.log10(mVaZa)
            if rat==float('inf'):
                final_ZV_ratio[i0,j0]=float('NaN')
            else:
                final_ZV_ratio[i0,j0]=rat
                mortality_ratios[i0,j0]=np.mean(mVa*100/(mZa+mVa+d))
        else:
            final_ZV_ratio[i0,j0]=float("NaN")
            mortality_ratios[i0,j0]=float("NaN")

        if surv in states_with_VandZ:
            mVaZa=np.mean(Va*Qz/(Za*Qv))
            rat=np.log10(mVaZa)
            if rat==float('inf'):
                final_ZV_ratio_count[i0,j0]=float('NaN')
            else:
                final_ZV_ratio_count[i0,j0]=rat
        else:
            final_ZV_ratio_count[i0,j0]=float("NaN")

        if surv_Z:
            final_Z[i0,j0] = np.log10(np.mean(Za))
        else:
            final_Z[i0,j0]= float("NaN")

        if surv_V:
            final_V[i0,j0] = np.log10(np.nanmean(Va))
        else:
            final_V[i0,j0] = float("NaN")

        if surv_S or surv_R or surv_V or surv_I:
            final_R[i0,j0] = np.log10(np.mean(Ra))
        else:
            final_R[i0,j0] = float("NaN")
        
        final_Surv[i0,j0]=surv

        # comparison to realistic and target concentrations
        mean_A=np.nanmean(Sa+Ra)/Qp
        mean_V=np.nanmean(Va)/Qv
        mean_Z=np.nanmean(Za)/Qz
        vals_r=[mean_A, mean_V, mean_Z]
        if surv in [5,6,7]:

            if mean_A > A_cond_low and mean_A < A_cond_high and mean_V>V_cond_low and mean_V<V_cond_high and mean_Z>Z_cond_low and mean_Z<Z_cond_high:
                valid_concentrations[i0,j0]=1
            if valid_concentrations[i0,j0]==1 and mortality_ratios[i0,j0]>perc_cond_low and mortality_ratios[i0,j0]<perc_cond_high:
                valid_concentrations_bis[i0,j0]=1
            if valid_concentrations[i0,j0]==1 and mean_V/mean_A>=1:
                valid_concentrations_ter[i0,j0]=1
            if valid_concentrations[i0,j0]==1 and valid_concentrations_bis[i0,j0]==1 and valid_concentrations_ter[i0,j0]==1:
                valid_concentrations_4[i0,j0]=1

            distance_to_t = absolute_error(target_conc, vals_r)

            distance_to_target[i0,j0]=distance_to_t
        else:
            distance_to_target[i0,j0]=float("NaN")

        if surv_V  or surv_S or surv_R:
          distance_to_target_A[i0,j0]=vals_r[0]-target_conc[0]
        else:
          distance_to_target_A[i0,j0]=float("NaN")
        if surv_V:
          distance_to_target_V[i0,j0]=vals_r[1]-target_conc[1]
        else:
          distance_to_target_V[i0,j0]=float("NaN")
        if surv_Z:
          distance_to_target_Z[i0,j0]=vals_r[2]-target_conc[2]
        else:
          distance_to_target_Z[i0,j0]=float("NaN")

    if spec_value!='no':
        aticky=[]
        alabel_ticky=[]
        atickx=[i*9 if i>0 else 0 for i in range(5) ]
        alabel_tickx= [phis[i] for i in atickx]
    else:
        atickx=[i*9 if i>0 else 0 for i in range(5) ]
        aticky=[i*9+9 if i!=0 else i*9+9 for i in range(3,-1,-1)]
        alabel_tickx= [phis[i] for i in atickx]
        atickyr=[i*9+8 if i!=0 else i*9+9 for i in range(0,4)]
        alabel_ticky= [phi_over_phir[i-9] for i in atickyr]
        aticky[0]=aticky[0]-1

    final_S[final_S==-inf]=float("NaN")
    final_Z[final_Z==-inf]=float("NaN")
    final_V[final_V==-inf]=float("NaN")
    final_I[final_I==-inf]=float("NaN")
    final_R[final_R==-inf]=float("NaN")
    final_perc_inf[final_perc_inf==-inf]=float("NaN")
    final_perc_res[final_perc_res==-inf]=float("NaN")
    final_pers[final_pers==-inf]=float("NaN")
    final_mods[final_mods==-inf]=float("NaN")
    
    # time series in the case of coexistence with dz2==0
    if dz2==0:
        pp = PdfPages('SVRZ_model_phi_latent_period_time_series_coex_'+suffix+'.pdf')
        peaks_to_save=[]
        if len(to_plot)>50:
            to_plot=to_plot[0::10]
        for k in to_plot:
            tit='phi='+str(k[0])+' phi/phir='+str(k[1])
            phi=k[0]
            if param=='phir':
                if k[1] !=0:
                    phir=phi/k[1]
                else:
                    phir=0
                epsr=eps
            elif param=='epsr':
                if k[1] !=0:
                    epsr=eps/k[1]
                else:
                    epsr=0
                phir=phi
            nyears0=20
            ndays=365*nyears0
            result=simulation_SVRZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur, dt, ndays, init_conditions)
            S=result[0]
            V=result[1]
            R=result[2]
            Z=result[3]
            i1=0
            i2=len(S)-1
            make_1_plot(S, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
            i1=0
            i2=round(365*5/dt)
            make_1_plot(S, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
            
            i1=round(365*(nyears0-1)/dt)
            i2=len(S)-1
            make_1_plot(S, V,R, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)

            S_sub=np.array(S[i1:i2])
            V_sub=np.array(V[i1:i2])
            R_sub=np.array(R[i1:i2])
            Z_sub=np.array(Z[i1:i2])

            S_peak, S_peak_l, S_peak_b=find_low_and_high(S_sub)
            V_peak, V_peak_l, V_peak_b=find_low_and_high(V_sub)
            R_peak, R_peak_l, R_peak_b=find_low_and_high(R_sub)
            Z_peak, Z_peak_l, Z_peak_b=find_low_and_high(Z_sub)
            peaks_to_save.append([phi, phir, epsr,S_peak*dt, S_peak_l*dt, V_peak*dt, V_peak_l*dt, R_peak*dt, R_peak_l*dt,Z_peak*dt, Z_peak_l*dt])
            if not np.isnan(S_peak) and not np.isnan(S_peak_b):
                make_1_plot(S_sub, V_sub, R_sub,Z_sub, int(S_peak), int(S_peak_b), tit,dt,pp, Qp, Qv, Qz)
        peaks_to_save=np.array(peaks_to_save)
        write_matrix(peaks_to_save, 'peaks_SVRZ_'+suffix+'.txt', sep=' ')
        pp.close()

    # time series at given adsorption rates for resistance strength == 10 (phi/phir)
    if indice in [1,2,3]:
        phis_to_plot=[1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    else:
        phis_to_plot=[1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
    pp = PdfPages(type_m+'_model_phi_latent_period_time_series_'+suffix+'.pdf')
    for phi in phis_to_plot:
        if param=='phir':
            phir=phi/10
            epsr=1
        else:
            phir=phi
            epsr=eps/10
        nyears0=5
        ndays=365*nyears0
        result=simulation_SVRZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur, dt, ndays, init_conditions)
        S=result[0]
        V=result[1]
        R=result[2]
        Z=result[3]

        i1=0
        i2=len(S)-1
        tit='phi='+str(phi)
        make_plots(S, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
            
        i3=round((365)/dt)
        make_plots(S, V, R,Z, i1, i3, tit,dt,pp, Qp, Qv, Qz)
    pp.close()

    # time series at given adsorption rates for fully resistant type
    pp = PdfPages(type_m+'_model_phi_latent_period_time_series_full-res_'+suffix+'.pdf')
    for phi in phis_to_plot:
        if param=='phir':
            phir=0
            epsr=1
        else:
            phir=phi
            epsr=0
        nyears0=5
        ndays=365*nyears0
        result=simulation_SVRZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur, dt, ndays, init_conditions)
        S=result[0]
        V=result[1]
        R=result[2]
        Z=result[3]

        i1=0
        i2=len(S)-1
        tit='phi='+str(phi)
        make_plots(S, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
        i3=round((365)/dt)
        make_plots(S, V, R,Z, i1, i3, tit,dt,pp, Qp, Qv, Qz)
    pp.close()

    # main pdf with matrices over the parameter space
    zeros_mat=np.zeros((len(phis),len(phi_over_phir)))
    pp = PdfPages(type_m+'_model_phi_versus_phir_'+suffix+'.pdf')
    if spec_value=='no':
        if param=='phir':
            ylb='phi/phir'
        else:
            ylb='eps/epsr'
    else:
        ylb=''
    
    mi=0
    mx=10
    bounds=[0,1,2,3,4,5,6,7,8,9,10]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=0
    mx=n_state-1
    colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx+2)])
    cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
    colors.append(cn)
    colors[10]=(0.8039, 0.5216, 0.2471)
    colors=tuple(colors)
    cmap = matplotlib.colors.ListedColormap(colors)
    axi=plot_with_scale(final_state_all,cmap ,mi,mx,atickx, aticky, alabel_tickx, alabel_ticky, 'State reached', norm=norm, yl=ylb)
    coex_mat=np.zeros( (len(phis),len(phi_over_phir)) ) 
    coex_mat[final_state_all>7]=1
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


    cmap=matplotlib.colormaps['bwr']
    rcmap = cmap.reversed()
    absmax=max([abs(np.nanmin(np.array(final_ZV_ratio))), abs(np.nanmax(np.array(final_ZV_ratio)))])
    mi=-absmax
    mx=absmax
    plot_with_scale_bis(final_ZV_ratio, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Grazer)', yl=ylb)
    pp.savefig()


    absmax=max([abs(np.nanmin(np.array(final_ZV_ratio_count))), abs(np.nanmax(np.array(final_ZV_ratio_count)))])
    mi=-absmax
    mx=absmax
    print(absmax)
    plot_with_scale_bis(final_ZV_ratio_count, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Grazer) count', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z))
    mx=np.nanmax(np.array(final_Z))
    plot_with_scale(final_Z, 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Zooplankton (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z)-mt.log10(Qz))
    mx=np.nanmax(np.array(final_Z)-mt.log10(Qz))
    plot_with_scale(np.array(final_Z)-mt.log10(Qz), 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V))
    mx=np.nanmax(np.array(final_V))
    plot_with_scale(final_V, 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Virus (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V)-mt.log10(Qv))
    mx=np.nanmax(np.array(final_V)-mt.log10(Qv))
    plot_with_scale(np.array(final_V)-mt.log10(Qv), 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus (log10(ind.L-1))', yl=ylb)
    pp.savefig()


    mi=np.nanmin(np.array(final_S))
    mx=np.nanmax(np.array(final_S))
    plot_with_scale(final_S, 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Susceptible (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_S)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_S)-mt.log10(Qp))
    plot_with_scale(np.array(final_S)-mt.log10(Qp), 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible (log10(ind.L-1))', yl=ylb)
    pp.savefig()


    mi=np.nanmin(np.array(final_R))
    mx=np.nanmax(np.array(final_R))
    plot_with_scale(final_R, 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Resistant (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_R)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_R)-mt.log10(Qp))
    plot_with_scale(np.array(final_R)-mt.log10(Qp), 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Resistant (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=0
    mx=100
    plot_with_scale_bis(final_perc_res, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of resistant cells' , yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_pers))
    mx=np.nanmax(np.array(final_pers))
    plot_with_scale(final_pers, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Dominant period (days)', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_mods))
    mx=np.nanmax(np.array(final_mods))
    plot_with_scale(final_mods, 'inferno', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'FFT modulo', yl=ylb)
    pp.savefig()

    mi=0
    mx=1
    cmap = matplotlib.colors.ListedColormap(['green', 'red'])
    plot_with_scale(final_stab, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Equilibrium type', yl=ylb)
    pp.savefig()

    mi=0
    mx=400
    distance_to_target_bis=copy.deepcopy(distance_to_target)
    distance_to_target_bis[distance_to_target_bis>400]=400
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


    to_plot=[valid_concentrations, valid_concentrations_bis, valid_concentrations_ter, valid_concentrations_4]
    titls=['Valid concentrations', 'Valid concentrations + % virus kill', 'Valid concentrations +  V/A ratio', 'Valid concentrations + % virus kill + V/A ratio']
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

    if spec_value=='no' and not np.isnan(ind_min):
        phi_index, rat_index = np.unravel_index(ind_min, distance_to_target.shape)
        min_d=distance_to_target[phi_index, rat_index]
        print(rat_index)
        phi_min=phis[phi_index]
        rat_min=phi_over_phir[rat_index]

        distance_to_target_val=copy.deepcopy(distance_to_target)
        distance_to_target_val[valid_concentrations_4==0]=np.nan

        try:
            ind_min_val=np.nanargmin(distance_to_target_val)
        except:
            ind_min_val=np.nan
            min_d_val=np.nan
            phi_min_val=np.nan
            rat_min_val=np.nan

        if not np.isnan(ind_min_val):
            phi_index_val, rat_index_val = np.unravel_index(ind_min_val, distance_to_target.shape)
            min_d_val=distance_to_target[phi_index_val, rat_index_val]
            phi_min_val=phis[phi_index_val]
            rat_min_val=phi_over_phir[rat_index_val]

            valid_concentrations_4[phi_index_val, rat_index_val]=2
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_index_par,:]=3
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            
            cond = (np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 10) & (np.nan_to_num(valid_concentrations_4, nan=np.inf) != 3)
            
            valid_concentrations_4[cond]=2
            vl=distance_to_target_val[valid_concentrations_4==3]<min_d_val+10
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_a_index_par]=4
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            mx=4
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +  % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

        else:
            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_index_par, :]=3
            valid_concentrations_4[phi_a_index_par, :]=4
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +% virus kill + V/A ratio', yl=ylb)
            pp.savefig()
        
        if 'vl' not in globals():
            vl=np.nan
        else:
            vl=vl[0]
        to_write=[min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val, vl]
        if nyears==20:
            write_vector(to_write, 'optimums_SVRZ_'+suffix+'.txt', ' ')
    elif not np.isnan(ind_min) and spec_value!='no':
        min_d=distance_to_target[ind_min]
        phi_min=phis[ind_min]

        try:
            distance_to_target_val=copy.deepcopy(distance_to_target)
            distance_to_target_val[valid_concentrations_4==0]=np.nan
            ind_min_val=np.nanargmin(distance_to_target_val)
        except:
            ind_min_val, min_d_val, phi_min_val=np.nan,np.nan,np.nan
        if not np.isnan(ind_min_val):
            phi_min_val=phis[ind_min_val]
            min_d_val=distance_to_target[ind_min_val]
            valid_concentrations_4[ind_min_val]=2

            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations  + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))

            valid_concentrations_4[phi_index_par]=3
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 10) & (np.nan_to_num(valid_concentrations_4, nan=np.inf) != 3)

            valid_concentrations_4[cond]=2
            
            vl=distance_to_target_val[valid_concentrations_4==3]<min_d_val+10
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations+ % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_a_index_par]=4
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            mx=4
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +  % virus kill + V/A ratio', yl=ylb)
            pp.savefig()
        else:
            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
            
            valid_concentrations_4[phi_index_par]=3
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_a_index_par]=4
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            mx=4
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +  % virus kill + V/A ratio', yl=ylb)
            pp.savefig()

        if 'vl' not in globals():
            vl=np.nan
        else:
            vl=vl[0]
        to_write=[min_d, phi_min, min_d_val, phi_min_val, vl]
        if nyears==20:
            write_vector(to_write, 'optimums_SVRZ_'+suffix+'.txt', ' ')
    elif np.isnan(ind_min):
        phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
        phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
        valid_concentrations_4[phi_index_par]=3
        valid_concentrations_4[phi_a_index_par]=4
        mi=0
        mx=np.nanmax(valid_concentrations_4)
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
        plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +% virus kill + V/A ratio', yl=ylb)
        pp.savefig()

        if spec_value=='no':
            min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val, vl=np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan
            to_write=[min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val, vl]
            if nyears==20:
                write_vector(to_write, 'optimums_SVRZ_'+suffix+'.txt', ' ')
        else:
            min_d, phi_min, min_d_val, phi_min_val, vl=np.nan, np.nan, np.nan, np.nan, np.nan
            to_write=[min_d, phi_min, min_d_val, phi_min_val, vl]
            if nyears==20:
                write_vector(to_write, 'optimums_SVRZ_'+suffix+'.txt', ' ')

    pp.close()

    if nyears==20:
        write_matrix(exclusion_times, 'exclusion_times_SVRZ_'+suffix+'.txt', ' ')
