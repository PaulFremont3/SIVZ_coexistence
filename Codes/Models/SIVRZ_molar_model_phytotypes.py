import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from numpy import inf
import copy
from SIVRZ_functions import *
from scipy.signal import find_peaks
from generic_functions import *


if __name__ == '__main__':
  indice=int(sys.argv[1]) # phytoplankton type, 0: small diatom, 1, picoeukaryote, 2, synechococcus, 3, prochlorococcus
  eff=sys.argv[2] # loss of infection rate (0)
  gr_ratio=sys.argv[3] # cost of resistance
  type_m=str(sys.argv[4]) # model type: SIVRZ (obsolete)
  otype=str(sys.argv[5]) # ocean type
  dz2=float(sys.argv[6]) # quadratic mortality term of zooplankton
  m2=float(sys.argv[7]) # quadratic mortality term of virus
  param=str(sys.argv[8]) # param to explore: lp_phir, phir, lp_epsr, epsr. If lp_epsr or lp_phir => latent period parameter space is explored, else: resistance strength 
  phi_ratio=float(sys.argv[9]) # specify fixed resistance strength in case param==lp_epsr or param==lp_phir. if 0 => full resistance

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
  # Z quota
  Qz=Q_grazer(rz)

  # V quota
  r_virus=[20,80,35,35]
  Qvs=[Q_virus(r) for r in r_virus]

  Qps=[]
  for i in range(100):
    Qps.append(Q_diatom(Vols[i]))
  for i in range(100, 200):
    Qps.append(Q_eukaryotes(Vols[i]))
  for i in range(200, 400):
    Qps.append(Q_cyanobacteria(Vols[i]))

  # choose the 4 smallest cells of each type (diatom, euk, syn and pro)
  N=4 # number of algae
  if N==4:
      ind1=np.argmin(np.array(Qps[0:100]))
      ind2=np.argmin(np.array(Qps[100:200]))
      ind3=np.argmin(np.array(Qps[200:300]))
      ind4=np.argmin(np.array(Qps[300:400]))

      # molar nitrogen quotas
      Qp1=Qps[0:100][ind1]
      Qp2=Qps[100:200][ind2]
      Qp3=Qps[200:300][ind3]
      Qp4=Qps[300:400][ind4]

      # burst size
      bs1=round(betas[ind1],-1) 
      bs2=round(betas[100:200][ind2],-1)
      bs3=round(betas[200:300][ind3],-1)
      bs4=round(betas[300:400][ind4],-1)
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

      # nurtien affinity constant
      Nc1=Ncs[0:100][ind1]
      Nc2=Ncs[100:200][ind2]
      Nc3=Ncs[200:300][ind3]
      Nc4=Ncs[300:400][ind4]
      Ncs=[Nc1, Nc2, Nc3, Nc4]

      # volumes
      V1=Vols[0:100][ind1]
      V2=Vols[100:200][ind2]
      V3=Vols[200:300][ind3]
      V4=Vols[300:400][ind4]
      Vs=[V1, V2, V3, V4]

  # choose phyto type
  typePhytos=['Diatom', 'Eukaryote', 'Synechococcus', 'Prochlorochoccus']
  typePhyto=typePhytos[indice]

  bs=betas[indice]
  lp=lpes[indice]
  lat_per=lpes[indice]

  print(mu_max)
  print(Ncs)

  
  # suffix to identify output pdfs
  suffix=typePhyto+'_LP'+str(lp)+'_BS'+str(bs)+'_LOI'+eff+'_GR-R'+gr_ratio
  if dz2==0:
      suffix+='_no-dz2'
  if m2 !=0:
      suffix+='_m2-'+str(m2)
  if param=='lp_phir':
      suffix+='_lp_phi-ratio-'+str(phi_ratio)
  elif param=='lp_epsr':
      suffix+='_lp_eps-ratio-'+str(phi_ratio)
  elif param=='epsr':
      suffix+='_epsr'
  elif param=='phir':
      suffix+='_phir'

  # temperature dependency
  Temp=20
  R=0.5
  T_dep=1
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
      elif otype=='mesotrophic':
        Temp=20
        T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
        R=1
        mu_max=[mu_max[i]*R/(R+Ncs[i])*T_dep for i in range(N)]
        print(mu_max)
      elif otype=='oligotrophic':
        Temp=25
        T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
        R=0.1
        mu_max=[mu_max[i]*R/(R+Ncs[i])*T_dep for i in range(N)]
        print(mu_max)

  # theoretical virus phi
  virus_radius=r_virus[indice]
  Host_volume=Vs[indice]
  phi_e, phi_e_a=phivs_encounter(virus_radius, Host_volume, Temp,indice)

  # simulation params
  dt=1/48
  nyears=20
  ndays=round(365*nyears)

  # model params
  mu=mu_max[indice]
  Nc=Ncs[indice]
  mui=0 # infected phyto max growth rate
  eta=float(lp) #latent period in days
  beta=float(bs) # burst size
  print('lp and bs:')
  print(lp)
  print(bs)
  ds=0.1*T_dep # phyto mortality rate
  m=0.1*T_dep #virus mortality
  m2=m2*T_dep #virus qudratic mortality
  Qv=Qvs[indice] #nitrogen cell quota virus umolN.virus-1
  Qp=Qps[indice]
  eps=1 # adsorption efficiency coefficient
  epso=float(eff) # loss of infectivity rate
  eps_r=1e-6 # conversion to resistant type
  eps_lr=1e-6 # loss of resitance
  

  mur=mu*float(gr_ratio)
  
  
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
  CC=KC_s/(mu-ds)

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

  # explored parameter space
  if indice in [1,2,3]:
    phis=list(np.logspace(-12, -8, (12-8)*9+1))
  elif indice in [0]:
    phis=list(np.logspace(-11, -7, (12-8)*9+1))
  if param=='phir' or param=='epsr':
    phi_over_phir=list(np.logspace(0, 4, 4*9+1))
    phi_over_phir=phi_over_phir[0:len(phi_over_phir)-1]
    phi_over_phir.append(0)
    del phi_over_phir[0]
    print(phi_over_phir)
    grid_test = [(ph, ratio) for ph in phis for ratio in phi_over_phir]
    gsize=len(phi_over_phir)
  elif param=='lp_phir' or param=='lp_epsr':
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]
    grid_test = [(ph, lp) for ph in phis for lp in lps]
    gsize=len(lps)

  # matrices to fill outputs (see SIVZ file for specs)
  limit_val_mat = np.zeros( (len(phis),gsize) )

  mortality_ratios=np.zeros( (len(phis),gsize)) 
  mortality_ratios_Z = np.zeros( (len(phis),gsize) )
  mortality_ratios_O = np.zeros( (len(phis),gsize) )
  mortality_ratios_th = np.zeros( (len(phis),gsize) )
  mortality_ratios_Z_th = np.zeros( (len(phis),gsize) )
  mortality_ratios_O_th = np.zeros( (len(phis),gsize) )
  mortality_ratios_th[:]=np.nan
  mortality_ratios_Z_th[:]=np.nan
  mortality_ratios_O_th[:]=np.nan

  final_state_all=np.empty( (len(phis),gsize))  
  final_state_all[:]=np.nan

  exclusion_times=np.zeros( (len(phis),gsize) )

  # realistic concentrations and targets
  A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc=concentration_ranges(indice, otype)

  valid_concentrations=np.zeros( (len(phis),gsize) )
  valid_concentrations_bis=np.zeros( (len(phis),gsize) )
  valid_concentrations_ter=np.zeros( (len(phis),gsize) )
  valid_concentrations_4=np.zeros( (len(phis),gsize) )
  valid_concentrations_5=np.zeros( (len(phis),gsize) )

  distance_to_target=np.zeros( (len(phis),gsize) )
  distance_to_target_A=np.zeros( (len(phis),gsize) )
  distance_to_target_V=np.zeros( (len(phis),gsize) )
  distance_to_target_Z=np.zeros( (len(phis),gsize) )

  # theoretical stability
  n_state=11
  a_stable_states_matrix = [np.empty( (len(phis),gsize) ) for i in range(n_state)]
  o_stable_states_matrix = [np.empty( (len(phis),gsize) ) for i in range(n_state)]
  f_stable_states_matrix = [np.empty( (len(phis),gsize) ) for i in range(n_state)]
  for l in range(n_state):
        a_stable_states_matrix[l][:]=np.nan
        o_stable_states_matrix[l][:]=np.nan
        f_stable_states_matrix[l][:]=np.nan

  #theory matrices
  final_SIR_th=np.zeros( (len(phis),gsize) )
  final_S_th=np.zeros( (len(phis),gsize) )
  final_I_th=np.zeros( (len(phis),gsize) )
  final_R_th=np.zeros( (len(phis),gsize) )
  final_V_th=np.zeros( (len(phis),gsize) )
  final_Z_th=np.zeros( (len(phis),gsize) )
  final_perc_inf_th=np.zeros( (len(phis),gsize) )
  final_perc_res_th=np.zeros( (len(phis),gsize) )
  final_perc_res_tot_th=np.zeros( (len(phis),gsize) )
  final_ZV_ratio_th=np.zeros( (len(phis),gsize) )
  final_SIR_th[:]=np.nan
  final_S_th[:]=np.nan
  final_I_th[:]=np.nan
  final_R_th[:]=np.nan
  final_V_th[:]=np.nan
  final_Z_th[:]=np.nan
  final_perc_inf_th[:]=np.nan
  final_perc_res_th[:]=np.nan
  final_perc_res_tot_th[:]=np.nan
  final_ZV_ratio_th[:]=np.nan

  npp_matrix= np.zeros( (len(phis),gsize) )

  alphas=[0 , 0.1, 0.5, 1]
  final_SIR= np.zeros( (len(phis),gsize) )
  final_S = np.zeros( (len(phis),gsize) )
  final_Surv = np.zeros( (len(phis),gsize) )
  final_Surv_alpha = [np.zeros( (len(phis),gsize) ) for i in range(len(alphas))]
  final_Surv_alpha0 = [np.zeros( (len(phis),gsize) ) for i in range(len(alphas))]
  final_Surv_alpha1 = [np.zeros( (len(phis),gsize) ) for i in range(len(alphas))]
  final_ZV_ratio = np.zeros( (len(phis),gsize) )
  final_ZV_ratio_count = np.zeros( (len(phis),gsize) )
  final_Z=np.zeros( (len(phis),gsize) )
  final_V=np.zeros( (len(phis),gsize) )
  final_I=np.zeros( (len(phis),gsize) )
  final_perc_inf= np.zeros( (len(phis),gsize) )
  final_perc_res= np.zeros( (len(phis),gsize) )
  final_perc_res_tot= np.zeros( (len(phis),gsize) )
  final_R=np.zeros( (len(phis),gsize) )
  final_pers=np.zeros( (len(phis),gsize) )
  final_mods=np.zeros( (len(phis),gsize) )
  final_stab=np.empty( (len(phis),gsize) )
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

  cases=np.zeros( (len(phis),gsize) ) 

  # main loop over parameter space
  init_conditions=[0.001, 0.001, 0.001, 0.001, 0.001]
  for k in grid_test:
    phi=k[0]
    if param=='phir':
      if k[1] !=0:
        phir=phi/k[1]
        epsr=eps
      else:
        phir=0
        epsr=eps
    
      i0=np.where(np.array(phis)==phi)[0][0]
      j0=np.where(np.array(phi_over_phir)==k[1])[0][0]
    elif param=='epsr':
      if k[1] !=0:
        phir=phi
        epsr=eps/k[1]
      else:
        phir=phi
        epsr=0

      i0=np.where(np.array(phis)==phi)[0][0]
      j0=np.where(np.array(phi_over_phir)==k[1])[0][0]
    elif param=='lp_phir':
      lp=k[1]
      eta=lp
      epsr=1
      if phi_ratio!=0:
        phir=phi/phi_ratio
      else:
        phir=0
      i0=np.where(np.array(phis)==phi)[0][0]
      j0=np.where(np.array(lps)==k[1])[0][0]
    elif param=='lp_epsr':
      lp=k[1]
      eta=lp
      phir=phi
      if phi_ratio!=0:
        epsr=eps/phi_ratio
      else:
        epsr=0
      i0=np.where(np.array(phis)==phi)[0][0]
      j0=np.where(np.array(lps)==k[1])[0][0]

    # simulation
    alpha=1
    result=simulation_SIVRZ_rk4(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC,alpha,dt, ndays, init_conditions)
  
    S=result[0]
    I=result[1]
    V=result[2]
    R=result[3]
    Z=result[4]
    mZ=result[5]
    mV=result[6]
    npp=result[7]

    # limit of theory
    if phi*KC_s/Qp<10*ds and dz2!=0 and m2==0:
        limit_val_mat[i0,j0]=1
    
    V = list(map(lambda x: 0 if x == float('inf') else x, V))
    I = list(map(lambda x: 0 if x == float('inf') else x, I))
    S = list(map(lambda x: 0 if x == float('inf') else x, S))
    R = list(map(lambda x: 0 if x == float('inf') else x, R))
    Z = list(map(lambda x: 0 if x == float('inf') else x, Z))

    # theoretical stability
    eta=lp
    if dz2!=0:
        if m2==0:
          a_stable,  o_stable, f_stable=check_stable_states_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC)
        else:
          a_stable,  o_stable, f_stable=check_stable_states_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC)
        for l in range(n_state):
            a_stable_states_matrix[l][i0,j0]=a_stable[l]
            o_stable_states_matrix[l][i0,j0]=o_stable[l]
            f_stable_states_matrix[l][i0,j0]=f_stable[l]

    # theory
    if type_m=='SIVRZ' and eff=='0' and dz2!=0 and m2==0:
      #1st eq
      alpha=1
      S_star, I_star, V_star, R_star, Z_star=equilibrium1_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
      if V_star>0 and Z_star>0 and R_star>0 and S_star>0 and I_star>0:
        theor_surv_phi.append(i0)
        theor_surv_phi_over_phir.append(gsize-j0-1)
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, Z_star)
        st=0
        for ei in eigs:
          ei_r=ei.real
          if ei_r>0:
            st=1
        final_stab[i0,j0]=st
        if st==1:
          final_unstab_phis.append(i0)
          final_unstab_phi_over_phirs.append(gsize-j0-1)

        if S_star*100/(S_star+R_star) <1:
          theor_surv_phi_all.append(i0)
          theor_surv_phi_over_phir_all.append(gsize-j0-1)
      for l, alpha in enumerate(alphas):
        if alpha !=0:
          S_star2, I_star2, V_star2, R_star2, Z_star2=equilibrium1_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
          if V_star2>Qv and Z_star2>Qz and R_star2>Qp and S_star2>Qp and I_star2>Qp:
            final_Surv_alpha[l][i0,j0]=1
        elif alpha==0:
          S_star2, I_star2, V_star2, R_star2, Z_star2=equilibrium1_SIVRZ_bis(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
          if V_star2>Qv and Z_star2>Qz and R_star2>Qp and S_star2>Qp and I_star2>Qp:
            final_Surv_alpha[l][i0,j0]=1


      #2nd eq
      alpha=1
      R_star0, I_star0, V_star0, S_star0, Z_star0=equilibrium2_SIVRZ(mur, mui, eta, beta, phir, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,phiz,eps_z,dz,dz2, CC, alpha)
      if V_star0>0 and Z_star0>0 and R_star0>0 and I_star0>0 and S_star0>0:
        theor_surv_phi0.append(i0)
        theor_surv_phi_over_phir0.append(gsize-j0-1)

      for l, alpha in enumerate(alphas):
        R_star0, I_star0, V_star0, S_star0, Z_star0=equilibrium2_SIVRZ(mur, mui, eta, beta, phir, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,phiz,eps_z,dz,dz2, CC, alpha)
        if V_star0>Qv and Z_star0>Qz and R_star0>Qp and I_star0>Qp and S_star0>0:
          final_Surv_alpha0[l][i0,j0]=1
            
      #3rd eq
      alpha=1
      S_star1, I_star1, V_star1, R_star1, Z_star1=equilibrium2_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
      if V_star1>0 and Z_star1>0 and S_star1>0 and I_star1>0 and R_star1>0:
        theor_surv_phi1.append(i0)
        theor_surv_phi_over_phir1.append(gsize-j0-1)

      for l, alpha in enumerate(alphas):
        S_star1, I_star1, V_star1, R_star1, Z_star1=equilibrium2_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)

        if V_star1>0 and Z_star1>0 and S_star1>0 and I_star1>0 and R_star1>0:
          final_Surv_alpha1[l][i0,j0]=1
    elif type_m=='SIVRZ' and eff=='0' and dz2!=0 and m2!=0:
      if eta<=7:
        limit_val_mat[i0,j0]=1
      #1st eq
      alpha=1
      S_star, I_star, V_star, R_star, Z_star,case =equilibrium1_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
      cases[i0,j0]=case-1
      if V_star>0 and Z_star>0 and S_star>0 and R_star>0  and I_star>0:
        theor_surv_phi.append(i0)
        theor_surv_phi_over_phir.append(gsize-j0-1)
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, Z_star)
        st=0
        for ei in eigs:
          ei_r=ei.real
          if ei_r>0:
            st=1
        final_stab[i0,j0]=st
        if st==1:
          final_unstab_phis.append(i0)
          final_unstab_phi_over_phirs.append(gsize-j0-1)

      #2nd eq (note the inversion of S and R parameters)
      alpha=1
      R_star0, I_star0, V_star0, S_star0, Z_star0=equilibrium2_SIVRZ_m2(mu=mur, mui=mui, lp=eta, beta=beta, phi=phir, d=ds, m=m,m2=m2, Qv=Qv, Qp=Qp,Qz=Qz,  eps=epsr,epsr=eps, epso=epso,eps_r=eps_r,eps_lr=eps_lr, phir=phi,mur=mu,phiz=phiz,eps_z=eps_z,dz=dz,dz2=dz2, CC=CC, alph=alpha)
      if V_star0>0 and Z_star0>0 and R_star0>0 and I_star0>0 and S_star0>0:
        theor_surv_phi0.append(i0)
        theor_surv_phi_over_phir0.append(gsize-j0-1)

      #3rd eq
      alpha=1
      S_star1, I_star1, V_star1, R_star1, Z_star1=equilibrium2_SIVRZ_m2(mu=mu, mui=mui, lp=eta, beta=beta, phi=phi, d=ds, m=m,m2=m2, Qv=Qv, Qp=Qp,Qz=Qz,  eps=eps, epsr=epsr,epso=epso,eps_r=eps_r,eps_lr=eps_lr,phir=phir,mur=mur,phiz=phiz,eps_z=eps_z,dz=dz,dz2=dz2, CC=CC, alph=alpha)
      if V_star1>0 and Z_star1>0 and S_star1>0 and I_star1>0 and R_star1>0:
        theor_surv_phi1.append(i0)
        theor_surv_phi_over_phir1.append(gsize-j0-1)    
    if type_m=='SIVRZ' and eff=='0' and dz2!=0:
      cond=V_star>0 and Z_star>0 and S_star>0 and R_star>0  and I_star>0
      cond0=V_star0>0 and Z_star0>0  and R_star0>0  and I_star0>0 and S_star0>0
      cond1=V_star1>0 and Z_star1>0 and S_star1>0 and R_star1>0  and I_star1>0

      if cond:
        S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star,R_star, I_star, V_star,Z_star
       
        fst=f_stable[0]
        fst0=f_stable[2]
        fst1=f_stable[1]
        
        if fst0 in [1,2] and fst==0:
          S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star0,R_star0, I_star0, V_star0,Z_star0
        elif fst1 in [1,2] and fst==0:
          S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star1,R_star1, I_star1, V_star1,Z_star1
    
      elif not cond and cond1 and not cond0:
        S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star1,R_star1, I_star1, V_star1,Z_star1
      elif not cond and cond0 and not cond1:
        S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star0,R_star0, I_star0, V_star0,Z_star0
      elif cond1 and cond0:
        fst0=f_stable[2]
        fst1=f_stable[1]

        if fst0 in [1,2] and fst1==0:
          S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star0,R_star0, I_star0, V_star0,Z_star0
        elif fst1 in [1,2] and fst0==0:
          S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star1,R_star1, I_star1, V_star1,Z_star1
        else:
          S_star_p, R_star_p, I_star_p, V_star_p,Z_star_p=S_star1,R_star1, I_star1, V_star1,Z_star1
       
      if cond0 or cond1 or cond:
        final_SIR_th[i0,j0]=mt.log10(S_star_p+I_star_p+R_star_p)
        final_S_th[i0,j0]=mt.log10(S_star_p)
        final_I_th[i0,j0]=mt.log10(I_star_p)
        final_R_th[i0,j0]=mt.log10(R_star_p)
        final_V_th[i0,j0]=mt.log10(V_star_p)
        final_Z_th[i0,j0]=mt.log10(Z_star_p)
        final_perc_inf_th[i0,j0]=I_star_p*100/(I_star_p+S_star_p+R_star_p)
        final_perc_res_th[i0,j0]=R_star_p*100/(R_star_p+S_star_p)
        final_perc_res_tot_th[i0,j0]=R_star_p*100/(R_star_p+S_star_p+I_star_p)
        final_ZV_ratio_th[i0,j0]=mt.log10(V_star_p/Z_star_p)

        mZt=phiz*Z_star_p
        pi=I_star_p/(I_star_p+S_star_p+R_star_p)
        mVt=pi/eta
        mortality_ratios_th[i0,j0]=mVt*100/(mZt+mVt+ds)
        mortality_ratios_Z_th[i0,j0]=mZt*100/(mZt+mVt+ds)
        mortality_ratios_O_th[i0,j0]=ds*100/(mZt+mVt+ds)
    

    # extract last year of simulation
    i1=int(365*(nyears-1)/dt)+1
    i2=len(V)
    Sa=np.array(S[i1:i2])
    Ia=np.array(I[i1:i2])
    Ra=np.array(R[i1:i2])
    Va=np.array(V[i1:i2])
    Za=np.array(Z[i1:i2])

    mZa=np.array(mZ[i1:i2])
    mVa=np.array(mV[i1:i2])

    NPPa=np.array(npp[i1:i2])

    # existence criterion
    surv_S=np.max(Sa/Qp)>1 #and np.log10(np.nanmean(Sa))> -5 #or np.nanmean(Sa/Qv)>1e-50
    surv_R=np.max(Ra/Qp)>1 #and np.log10(np.nanmean(Ra))> -5
    surv_V = np.max(Va/Qv)>1 #and np.log10(np.nanmean(Va/Qv))>4 #or np.nanmean(Va/Qv) > 1e-50
    surv_Z=np.max(Za/Qz)>1 #and np.log10(np.nanmean(Za/Qz))>1# #or np.nanmean(Za/Qz) > 1e-50
    surv_I=np.max(Ia/Qp)>1 #and np.log10(np.nanmean(Ia/Qp))>2

    # fill matrices
    if surv_V or surv_I or surv_S or surv_R:
      final_S[i0,j0]=np.log10(np.mean(Sa))
      final_SIR[i0,j0]=np.log10(np.mean(Sa+Ia+Ra))
      final_perc_res[i0,j0]= np.mean(Ra*100/(Ra+Sa))
      final_perc_res_tot[i0,j0]= np.mean(Ra*100/(Ra+Sa+Ia))
      npp_matrix[i0,j0]=np.log10(np.mean(NPPa*106/16))
    else:
      final_S[i0,j0]=float("NaN")
      final_SIR[i0,j0]=float("NaN")
      final_perc_res[i0,j0]=float("NaN")
      final_perc_res_tot[i0,j0]=float("NaN")
      npp_matrix[i0,j0]=float("NaN")
  
    surv_S_bis = np.mean(Sa) > 0.00001 #or np.max(Sa/Qp)>1
    surv_R_bis = np.mean(Ra) > 0.00001  #or np.max(Ra/Qp)>1
    if surv_V:
      f0, modulos=fft_data(Va,dt)
      mx_m=np.max(modulos)
      mx_per=1/f0[np.nanargmax(modulos)]
      final_mods[i0,j0]=mt.log10(mx_m)
      final_pers[i0,j0]=mt.log10(mx_per)
    else:
      final_mods[i0,j0]=float("NaN")
      final_pers[i0,j0]=float("NaN")
      
    # final state reached by the system
    surv=np.nan
    if not surv_V  and not surv_Z and not surv_I and not surv_S and not surv_R:
      final_state_all[i0,j0]=0
      surv=0
    elif surv_R and not surv_V  and not surv_Z and not surv_I and not surv_S_bis:
      final_state_all[i0,j0]=1
      surv=1
    elif surv_S and not surv_R_bis and not surv_V  and not surv_Z and not surv_I:
      final_state_all[i0,j0]=2
      surv=1
    elif surv_R and not (surv_V or surv_I)  and surv_Z and not surv_S_bis:
      final_state_all[i0,j0]=3
      surv=2
    elif surv_S and not surv_R_bis and not (surv_V or surv_I)  and surv_Z:
      final_state_all[i0,j0]=4
      surv=2
    elif surv_R and  (surv_V or surv_I)  and not surv_Z and not surv_S_bis:
      final_state_all[i0,j0]=5
      surv=3
    elif surv_S and  (surv_V or surv_I)  and not surv_Z and not surv_R_bis:
      final_state_all[i0,j0]=6
      surv=3
    elif surv_S and  (surv_V or surv_I)  and not surv_Z and surv_R:
      final_state_all[i0,j0]=7
      surv=4
    elif not surv_S_bis and  (surv_V or surv_I)  and surv_Z and surv_R:
      final_state_all[i0,j0]=8
      surv=5
    elif surv_S and  (surv_V or surv_I)  and surv_Z and not surv_R_bis:
      final_state_all[i0,j0]=9
      surv=6
    elif surv_S and  (surv_V or surv_I)  and surv_Z and surv_R:
      final_state_all[i0,j0]=10
      surv=7

    if np.isnan(surv):
        surv=0
        final_state_all[i0,j0]=0
        print(surv_S)
        print(surv_V)
        print(surv_I)
        print(surv_Z)
        print(surv_R)
        print(surv_R_bis)
        print(surv_S_bis)

    if final_state_all[i0,j0]==10 and dz2==0:
      to_plot.append(k)
      print('here:')
      print(k)
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
    else:
      final_ZV_ratio[i0,j0]=float("NaN")
      
    if (surv_V or surv_I) and surv_Z:
      mortality_ratios[i0,j0]=np.mean(mVa*100/(mZa+mVa+ds))
      mortality_ratios_Z[i0,j0]=np.mean(mZa*100/(mZa+mVa+ds))
      mortality_ratios_O[i0,j0]=np.mean(ds*100/(mZa+mVa+ds))
    elif (surv_V or surv_I) and not surv_Z:
      mortality_ratios[i0,j0]=np.mean(mVa*100/(mVa+ds))
      mortality_ratios_O[i0,j0]=np.mean(ds*100/(mVa+ds))
      mortality_ratios_Z[i0,j0]=float("NaN")
    elif not (surv_V or surv_I) and surv_Z:
      mortality_ratios[i0,j0]=float("NaN")
      mortality_ratios_O[i0,j0]=np.mean(ds*100/(mZa+ds))
      mortality_ratios_Z[i0,j0]=np.mean(mZa*100/(mZa+ds))
    else:
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
      final_Z[i0,j0] = np.log10(np.nanmean(Za))
    else:
      final_Z[i0,j0]= float("NaN")

    if surv_V or surv_I:
      final_V[i0,j0] = np.log10(np.nanmean(Va))
    else:
      final_V[i0,j0] = float("NaN")

    if surv_V or surv_I:
      final_I[i0,j0] = np.log10(np.mean(Ia))
    else:
      final_I[i0,j0] = float("NaN")
    if surv_S or surv_R or surv_V or surv_I:
      final_R[i0,j0] = np.log10(np.mean(Ra))
    else:
      final_R[i0,j0] = float("NaN")
    if surv_V or surv_I:
      final_perc_inf[i0,j0]= np.mean(Ia*100/(Ia+Ra+Sa))
    else:
      final_perc_inf[i0,j0]=float("NaN")
    final_Surv[i0,j0]=surv

    # comparison to realistic concentration and ddistance to target concentrations
    mean_A=np.nanmean(Sa+Ia+Ra)/Qp
    mean_V=np.nanmean(Va)/Qv
    mean_Z=np.nanmean(Za)/Qz
    vals_r=[mean_A, mean_V, mean_Z]
    if surv in [5,6,7]:
      if mean_A > A_cond_low and mean_A < A_cond_high and mean_V>V_cond_low and mean_V<V_cond_high and mean_Z>Z_cond_low and mean_Z<Z_cond_high:
        valid_concentrations[i0,j0]=1
      if valid_concentrations[i0,j0]==1 and final_perc_inf[i0,j0]>I_cond_low and final_perc_inf[i0,j0]<I_cond_high:
        valid_concentrations_bis[i0,j0]=1
      if valid_concentrations[i0,j0]==1 and mortality_ratios[i0,j0]>perc_cond_low and mortality_ratios[i0,j0]<perc_cond_high:
        valid_concentrations_ter[i0,j0]=1
      if valid_concentrations[i0,j0]==1 and mean_V/mean_A>=1:
        valid_concentrations_4[i0,j0]=1
      if valid_concentrations_bis[i0,j0]==1 and valid_concentrations_ter[i0,j0]==1 and valid_concentrations_4[i0,j0]==1:
        valid_concentrations_5[i0,j0]=1
      distance_to_t =  absolute_error(target_conc, vals_r)
    
      distance_to_target[i0,j0]=distance_to_t
    else:
      distance_to_target[i0,j0]=float("NaN")

    if surv_V or surv_I or surv_R or surv_S:
      distance_to_target_A[i0,j0]=vals_r[0]-target_conc[0]
    else:
      distance_to_target_A[i0,j0]=float("NaN")
    if surv_V or surv_I:
      distance_to_target_V[i0,j0]=vals_r[1]-target_conc[1]
    else:
      distance_to_target_V[i0,j0]=float("NaN")
    if surv_Z:
      distance_to_target_Z[i0,j0]=vals_r[2]-target_conc[2]
    else:
      distance_to_target_Z[i0,j0]=float("NaN")

  if param=='phir' or param=='epsr':
    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+9 if i!=0 else i*9+9 for i in range(3,-1,-1)]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=[i*9+8 if i!=0 else i*9+9 for i in range(0,4)]
    alabel_ticky= [phi_over_phir[i-9] for i in atickyr]
    aticky[0]=aticky[0]-1
    ylb='phi/phir'
    if param=='epsr':
        ylb='eps/epsr'
  elif param=='lp_phir' or param=='lp_epsr':
    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]
    ylb='latent period'

  final_S[final_S==-inf]=float("NaN")
  final_Z[final_Z==-inf]=float("NaN")
  final_V[final_V==-inf]=float("NaN")
  final_I[final_I==-inf]=float("NaN")
  final_R[final_R==-inf]=float("NaN")
  final_perc_inf[final_perc_inf==-inf]=float("NaN")
  final_perc_res[final_perc_res==-inf]=float("NaN")
  final_perc_res_tot[final_perc_res_tot==-inf]=float("NaN")
  final_pers[final_pers==-inf]=float("NaN")
  final_mods[final_mods==-inf]=float("NaN")

  
  print(to_plot)
  # time series in case of coexistence
  if dz2==0 and (param=='phir' or param=='epsr'):
    epsr=eps
    pp = PdfPages(type_m+'_model_phi_latent_period_time_series_coex_'+suffix+'.pdf')
    for k in to_plot[0::20]:
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
      result=simulation_SIVRZ_rk4(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC,alpha,dt, ndays, init_conditions)
      S=result[0]
      I=result[1]
      V=result[2]
      R=result[3]
      Z=result[4]
      i1=0
      i2=len(S)-1
      make_1_plot_SIVRZ(S, I, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
      i1=0
      i2=round(365*5/dt)
      make_1_plot_SIVRZ(S, I, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
    pp.close()

  # time series for specific values of adsorption rate with resistance strength =10 (phi/phir)
  if param=='phir' or param=='epsr':
    epsr=eps
    lat_per=lpes[indice]
    if indice in [1,2,3]:
      phis_to_plot=[1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    else:
      phis_to_plot=[1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
    pp = PdfPages(type_m+'_model_phi_latent_period_time_series_'+suffix+'.pdf')
    for phi in phis_to_plot:
      if  param=='phir': 
        phir=phi/10
        epsr=eps
      elif param=='epsr':
        phir=phi
        epsr=eps/10
      nyears0=5
      ndays=365*nyears0
      result=simulation_SIVRZ_rk4(mu, mui, lat_per, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC,alpha,dt, ndays, init_conditions)
      S=result[0]
      I=result[1]
      V=result[2]
      R=result[3]
      Z=result[4]
  
  
      i1=0
      i2=len(S)-1
      tit='phi='+str(phi)
      make_plots_SIVRZ(S, I, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
      i3=round((365)/dt)
      make_plots_SIVRZ(S, I, V, R,Z, i1, i3, tit,dt,pp, Qp, Qv, Qz)
    pp.close()

    # time series for specific values of adsorption rate with full resistance phir=0
    pp = PdfPages(type_m+'_model_phi_latent_period_time_series_full-res_'+suffix+'.pdf')
    for phi in phis_to_plot:
      if param=='phir':
        phir=0
        epsr=eps
      else:
        epsr=0
        phir=phi
      nyears0=5
      ndays=365*nyears0
      result=simulation_SIVRZ_rk4(mu, mui, lat_per, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC,alpha,dt, ndays, init_conditions)
      S=result[0]
      I=result[1]
      V=result[2]
      R=result[3]
      Z=result[4]


      i1=0
      i2=len(S)-1
      tit='phi='+str(phi)
      make_plots_SIVRZ(S, I, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
      i3=round((365)/dt)
      make_plots_SIVRZ(S, I, V, R,Z, i1, i3, tit,dt,pp, Qp, Qv, Qz)
    pp.close()


  # main pdf with output matrices over the parameter space
  zeros_mat=np.zeros((len(phis),gsize))
  if param=='phir':
    pp = PdfPages(type_m+'_model_phi_versus_phir_'+suffix+'.pdf')
  elif param=='epsr':
    pp = PdfPages(type_m+'_model_eps_versus_epsr_'+suffix+'.pdf')    
  elif param=='lp_phir' or param=='lp_epsr':
    pp = PdfPages(type_m+'_model_phi_versus_latent_period_'+suffix+'.pdf')

  bounds=[0,1,2,3,4,5,6,7,8,9,10]
  norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)
  mi=0
  mx=n_state-1
  colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx+2)])
  cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
  colors.append(cn)
  colors[10]=(0.8039, 0.5216, 0.2471)
  colors=tuple(colors)
  cmap = matplotlib.colors.ListedColormap(colors)

  axi=plot_with_scale(final_state_all,cmap ,mi,mx,atickx, aticky, alabel_tickx, alabel_ticky, 'State reached', norm=norm, yl=ylb)
  coex_mat=np.zeros( (len(phis),gsize) )
  coex_mat[final_state_all>7]=1
  coex_mat=np.transpose(coex_mat)
  coex_mat=np.flipud(coex_mat)
  draw_borders(coex_mat, 1.5, 0.05, 'white', axi)
  pp.savefig()


  if param=='lp_phir' or param=='lp_epsr':
    phi_index_par,phi_a_index_par, lp_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(lat_per-np.array(lps)))
    final_state_all_mod=copy.deepcopy(final_state_all)
    final_state_all_mod[phi_index_par, lp_index_par]=11
  else:
    phi_index_par,phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
    final_state_all_mod=copy.deepcopy(final_state_all)
    final_state_all_mod[phi_index_par, :]=11

  colors_bis=colors+('limegreen', 'cyan')
  cmap_R_bis= matplotlib.colors.ListedColormap(colors_bis)
  bounds_bis=[0,1,2,3,4,5,6,7,8,9,10,11,12]
  norm_bis = matplotlib.colors.BoundaryNorm(bounds_bis, cmap_R_bis.N-1)
  
  plot_with_scale(final_state_all_mod,cmap_R_bis,mi,mx+2, atickx, aticky, alabel_tickx, alabel_ticky, 'State reached + LHT model', norm=norm_bis, yl=ylb)
  pp.savefig()

  if param=='lp_phir' or param=='lp_epsr':
    final_state_all_mod[phi_a_index_par, lp_index_par]=12
  else:
    final_state_all_mod[phi_a_index_par,:]=12
  plot_with_scale(final_state_all_mod,cmap_R_bis,mi,mx+2, atickx, aticky, alabel_tickx, alabel_ticky, 'State reached + LHT models', norm=norm_bis, yl=ylb)
  pp.savefig()

  # theory
  if type_m=='SIVRZ' and eff=='0' and dz2!=0:
    mi=0
    mx=2
    cmap = matplotlib.colors.ListedColormap(['red', 'green','dodgerblue'])
    eq_types=['S,I,V,R,Z', 'S,I,V,0,Z', '0,I,V,R,Z', 'S,I,V,R,0', 'S,I,V,0,0', '0,I,V,R,0', 'S,0,0,0,Z', '0,0,0,R,Z', 'S,0,0,0,0', '0,0,0,R,0', '0,0,0,0,0']
    for l in range(n_state):
        matri=copy.deepcopy(f_stable_states_matrix[l])
        matri[limit_val_mat==0]=np.nan
        plot_with_scale_bis(matri, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, eq_types[l], yl=ylb)
        pp.savefig()

    # theory
    predicted_state=np.empty((len(phis), gsize))
    predicted_state[:]=np.nan
    jo=n_state-1
    to_check=list(range(10))
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
    
    alter_state=np.zeros((len(phis), gsize))
    u_state=np.zeros((len(phis), gsize))
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

    alter_state[limit_val_mat==0]=0
    u_state[limit_val_mat==0]=0

    alter_state[alter_state==1]=0
    alter_state[alter_state>=2]=1

    u_state[u_state!=n_state]=0
    u_state[u_state==n_state]=1

    mi=0
    mx=n_state-1
    predicted_state[limit_val_mat==0]=np.nan
    colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx+2)])
    cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
    colors.append(cn)
    colors[10]=(0.8039, 0.5216, 0.2471)
    colors=tuple(colors)
    cmap = matplotlib.colors.ListedColormap(colors)
    axi=plot_with_scale_bis(predicted_state,cmap ,mi,mx,atickx, aticky, alabel_tickx, alabel_ticky, 'Predicted state', yl=ylb)
    plt.contour(np.transpose(alter_state), colors=['limegreen', 'white'], levels=[0,1] , extent=[-1, len(phis), -1, gsize], origin='upper', linewidths=0.5)
    alter_state_bis=copy.deepcopy(alter_state)
    alter_state_bis[alter_state==1]=0
    alter_state_bis[alter_state==0]=1
    plt.contourf(np.transpose(alter_state_bis), colors=['limegreen', 'white'], levels=[0,0.9] , alpha=0.2, extent=[-1, len(phis), -1, gsize], origin='upper')

    plt.contour(np.transpose(u_state), colors=['black', 'white'], levels=[0,1] , extent=[-1, len(phis), -1, gsize], origin='upper', linewidths=0.5)
    u_state_bis=copy.deepcopy(u_state)
    u_state_bis[u_state==1]=0
    u_state_bis[u_state==0]=1
    plt.contourf(np.transpose(u_state_bis), colors=['black', 'white'], levels=[0,0.9] , alpha=0.2, extent=[-1, len(phis), -1, gsize], origin='upper')
    
    coex_mat=np.zeros( (len(phis),gsize) )
    coex_mat[predicted_state>7]=1
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
  plot_with_scale_bis(final_ZV_ratio_count, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Grazer) count', yl=ylb)
  pp.savefig()

  final_ZV_ratio_th[limit_val_mat==0]=np.nan
  absmax=max([abs(np.nanmin(np.array(final_ZV_ratio_th))), abs(np.nanmax(np.array(final_ZV_ratio_th)))])
  mi=-absmax
  mx=absmax
  plot_with_scale_bis(final_ZV_ratio_th, rcmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'log10(Virus/Grazer) th', yl=ylb)
  pp.savefig()
  
  mi=np.nanmin(np.array(final_Z))
  mx=np.nanmax(np.array(final_Z))
  plot_with_scale(final_Z, 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Zooplankton (log10(umol.L-1))', yl=ylb)
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
  plot_with_scale(final_V, 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Virus (log10(umol.L-1))', yl=ylb)
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

  mi=np.nanmin(np.array(final_SIR))
  mx=np.nanmax(np.array(final_SIR))
  plot_with_scale(final_SIR, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'S+I+R (log10(umol.L-1))', yl=ylb)
  pp.savefig()

  mi=np.nanmin(np.array(final_SIR)-mt.log10(Qp))
  mx=np.nanmax(np.array(final_SIR)-mt.log10(Qp))
  plot_with_scale(np.array(final_SIR)-mt.log10(Qp), 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I+R (log10(ind.L-1))', yl=ylb)
  pp.savefig()

  final_SIR_th[limit_val_mat==0]=np.nan
  mi=np.nanmin(np.array(final_SIR_th)-mt.log10(Qp))
  mx=np.nanmax(np.array(final_SIR_th)-mt.log10(Qp))
  plot_with_scale(np.array(final_SIR_th)-mt.log10(Qp), 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I+R th (log10(ind.L-1))', yl=ylb)
  pp.savefig()

  mi=np.nanmin(np.array(final_S))
  mx=np.nanmax(np.array(final_S))
  plot_with_scale(final_S, 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Susceptible (log10(umol.L-1))', yl=ylb)
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

  mi=np.nanmin(np.array(final_R))
  mx=np.nanmax(np.array(final_R))
  plot_with_scale(final_R, 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Resistant (log10(umol.L-1))', yl=ylb)
  pp.savefig()
  
  mi=np.nanmin(np.array(final_R)-mt.log10(Qp))
  mx=np.nanmax(np.array(final_R)-mt.log10(Qp))
  plot_with_scale(np.array(final_R)-mt.log10(Qp), 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Resistant (log10(ind.L-1))', yl=ylb)
  pp.savefig()

  final_R_th[limit_val_mat==0]=np.nan
  mi=np.nanmin(np.array(final_R_th)-mt.log10(Qp))
  mx=np.nanmax(np.array(final_R_th)-mt.log10(Qp))
  plot_with_scale(np.array(final_R_th)-mt.log10(Qp), 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Resistant th (log10(ind.L-1))', yl=ylb)
  pp.savefig()
  
  mi=np.nanmin(np.array(npp_matrix))
  mx=np.nanmax(np.array(npp_matrix))
  plot_with_scale(np.array(npp_matrix), 'YlGn', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'NPP (log10(umolC.L-1.d-1))', yl=ylb)
  pp.savefig()

  mi=0
  mx=100
  plot_with_scale_bis(final_perc_inf, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of infected cells' , yl=ylb)
  pp.savefig()

  final_perc_inf_th[limit_val_mat==0]=np.nan
  plot_with_scale_bis(final_perc_inf_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of infected cells th' , yl=ylb)
  pp.savefig()

  plot_with_scale_bis(final_perc_res, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of resistant cells' , yl=ylb)
  pp.savefig()

  final_perc_res_th[limit_val_mat==0]=np.nan
  plot_with_scale_bis(final_perc_res_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of resistant cells th' , yl=ylb)
  pp.savefig()

  plot_with_scale_bis(final_perc_res_tot, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of resistant cells (tot)' , yl=ylb)
  pp.savefig()

  final_perc_res_th[limit_val_mat==0]=np.nan
  plot_with_scale_bis(final_perc_res_tot_th, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of resistant cells th (tot)' , yl=ylb)
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
  titls=['Valid concentrations', 'Valid concentrations+ % infected', 'Valid concentrations + % virus kill','Valid concentrations + V/A ratio', 'Valid concentrations + % infected + % virus kill + V/A ratio']
  mi=0
  mx=1
  cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue'])
  for k,mato in enumerate(to_plot):
    plot_with_scale(mato, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, titls[k], yl=ylb)
    pp.savefig()
  
  
  ind_min=np.nanargmin(distance_to_target)
  if not np.isnan(ind_min):
    if param in ['lp_epsr', 'lp_phir']:
      phi_index, lp_index = np.unravel_index(ind_min, distance_to_target.shape)
      min_d=distance_to_target[phi_index, lp_index]
      phi_min=phis[phi_index]
      lp_min=lps[lp_index]

      distance_to_target_val=copy.deepcopy(distance_to_target)
      distance_to_target_val[valid_concentrations_5==0]=np.nan

      try:
        ind_min_val=np.nanargmin(distance_to_target_val)
      except:
        min_d_val=np.nan
        phi_min_val=np.nan
        lp_min_val=np.nan
        ind_min_val=np.nan

      if not np.isnan(ind_min_val):
        phi_index_val, lp_index_val = np.unravel_index(ind_min_val, distance_to_target.shape)
        min_d_val=distance_to_target_val[phi_index_val, lp_index_val]
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
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio')
        pp.savefig()
      
        if phi_a_index_par==phi_index_par:
          valid_concentrations_5[phi_a_index_par, lp_index_par]=3
        cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) > min_d_val + 10) & (np.nan_to_num(valid_concentrations_5, nan=np.inf) != 3) & (np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 50)
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan', 'darkorange'])
        valid_concentrations_5[cond]=5
        mx=5
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio')
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

      if 'vl' not in globals():
        vl=np.nan
      else:
        vl=vl[0]
      to_write=[min_d, phi_min,lp_min, min_d_val, phi_min_val, lp_min_val, vl]
      if nyears==20:
        write_vector(to_write, 'optimums_'+type_m+'_'+suffix+'.txt', ' ')
    elif param in ['epsr', 'phir']:
      phi_index, rat_index = np.unravel_index(ind_min, distance_to_target.shape)
      min_d=distance_to_target[phi_index, rat_index]
      phi_min=phis[phi_index]
      rat_min=phi_over_phir[rat_index]

      distance_to_target_val=copy.deepcopy(distance_to_target)
      distance_to_target_val[valid_concentrations_5==0]=np.nan

      try:
        ind_min_val=np.nanargmin(distance_to_target_val)
      except:
        min_d_val=np.nan
        phi_min_val=np.nan
        rat_min_val=np.nan
        ind_min_val=np.nan

      if not np.isnan(ind_min_val):
        phi_index_val, rat_index_val = np.unravel_index(ind_min_val, distance_to_target.shape)
        min_d_val=distance_to_target[phi_index_val, rat_index_val]
        phi_min_val=phis[phi_index_val]
        rat_min_val=phi_over_phir[rat_index_val]

        valid_concentrations_5[phi_index_val, rat_index_val]=2
        mi=0
        mx=2
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold'])
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
        pp.savefig()

        phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))

        valid_concentrations_5[phi_index_par]=3
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
        valid_concentrations_5[phi_a_index_par]=4
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
        mx=4
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio')
        pp.savefig()
        
        if phi_a_index_par==phi_index_par:
          valid_concentrations_5[phi_a_index_par]=3
        cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) > min_d_val + 10) & (np.nan_to_num(valid_concentrations_5, nan=np.inf) != 3) & (np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 50)
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan', 'darkorange'])
        valid_concentrations_5[cond]=5
        mx=5
        plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio')
        pp.savefig()

      if 'vl' not in globals():
        vl=np.nan
      else:
        vl=vl[0]
      to_write=[min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val]
      if nyears==20:
        write_vector(to_write, 'optimums_'+type_m+'_'+suffix+'.txt', ' ')
  else:
    if param in ['epsr', 'phir']:
      phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
      phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
      valid_concentrations_5[phi_index_par,:]=3
      valid_concentrations_5[phi_a_index_par,:]=4
      mi=0
      mx=4
      cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
      plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
      pp.savefig()
      
      min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val, vl=np.nan,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      to_write=[min_d, phi_min,rat_min, min_d_val, phi_min_val, rat_min_val, vl]
      if nyears==20:
        write_vector(to_write, 'optimums_'+type_m+'_'+suffix+'.txt', ' ')
    elif param in ['lp_epsr', 'lp_phir']:
      phi_index_par, lp_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis)))), np.nanargmin(np.absolute(lat_per-np.array(lps)))
      phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
      valid_concentrations_5[phi_index_par, lp_index_par]=3
      valid_concentrations_5[phi_a_index_par, lp_index_par]=3

      mi=0
      mx=4
      cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
      plot_with_scale(valid_concentrations_5, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % infected + % virus kill + V/A ratio', yl=ylb)
      pp.savefig()

      min_d, phi_min,lp_min, min_d_val, phi_min_val, lp_min_val, vl=np.nan,np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
      to_write=[min_d, phi_min,lp_min, min_d_val, phi_min_val, lp_min_val, vl]
      if nyears==20:
        write_vector(to_write, 'optimums_'+type_m+'_'+suffix+'.txt', ' ')

  pp.close()

  # save only if full simulation
  if nyears==20:
    write_matrix(exclusion_times, 'exclusion_times_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_S-mt.log10(Qp), 'model_data/Susceptible_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_I-mt.log10(Qp), 'model_data/Infected_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_R-mt.log10(Qp), 'model_data/Resistant_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_SIR-mt.log10(Qp), 'model_data/Tot_algae_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_V-mt.log10(Qv), 'model_data/Virus_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_Z-mt.log10(Qz), 'model_data/Zoop_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_ZV_ratio, 'model_data/ZV_ratio_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(distance_to_target, 'model_data/distances_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(distance_to_target_A, 'model_data/distance_A_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(distance_to_target_V, 'model_data/distance_V_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(distance_to_target_Z, 'model_data/distance_Z_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(valid_concentrations_5, 'model_data/validation_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_pers, 'model_data/Period_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(final_mods, 'model_data/FFT_modulo_'+type_m+'_'+suffix+'.txt', ' ')
    write_matrix(npp_matrix, 'model_data/NPP_'+type_m+'_'+suffix+'.txt', ' ')

  # save theoretical matrices
  write_matrix(final_S_th-mt.log10(Qp), 'model_data/Susceptible_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_I_th-mt.log10(Qp), 'model_data/Infected_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_SIR_th-mt.log10(Qp), 'model_data/Tot_algae_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_V_th-mt.log10(Qv), 'model_data/Virus_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_Z_th-mt.log10(Qz), 'model_data/Zoop_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_R_th-mt.log10(Qp), 'model_data/Resistant_th_'+type_m+'_'+suffix+'.txt', ' ')
  write_matrix(final_ZV_ratio_th, 'model_data/ZV_ratio_th_'+type_m+'_'+suffix+'.txt', ' ')
