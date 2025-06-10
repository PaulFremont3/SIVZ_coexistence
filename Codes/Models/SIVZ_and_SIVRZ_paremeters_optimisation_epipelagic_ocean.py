import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from numpy import inf
import copy
from SIVZ_functions import *
from SIVRZ_functions import *
from scipy.signal import find_peaks
from generic_functions import *
import os
import time
import psutil
import gc

if __name__ == '__main__':
    indice=int(sys.argv[1]) # phytoplankton type
    type_m=str(sys.argv[2]) # model type SIVZ, SIVZ_intra, SIVRZ or SIVRZ_intra
    POM=str(sys.argv[3]) # 0 or 1 
    chunk=str(sys.argv[4]) # chunk

    # lost life history traits
    Vols=load_vector('../trait_data/Vs_5.txt', sep=' ')
    Ncs=load_vector('../trait_data/Nc_dutkiewicz_5.txt', sep=' ')
    betas=load_vector('../trait_data/model_burst_size_nn-gam.txt', sep=' ')
    lps=load_vector('../trait_data/model_latent_period_nn-gam.txt', sep=' ')
    mu_max=load_vector('../trait_data/mumax_dutkiewicz_5.txt', sep=' ')
  
    if indice in [0,1]:
        rz=5
    elif indice in [2,3]:
        rz=2.5
    # zooplankton quota 
    Qz=Q_grazer(rz)

    # virus quota
    r_virus=[20,80,35,35]
    Qvs=[Q_virus(r) for r in r_virus]

    Qps=[]
    for i in range(100):
        Qps.append(Q_diatom(Vols[i]))
    for i in range(100, 200):
        Qps.append(Q_eukaryotes(Vols[i]))
    for i in range(200, 400):
        Qps.append(Q_cyanobacteria(Vols[i]))

    # choose the 4 smallest phytoplankton type of each subgroup (pro, syn, euk, diatom)
    N=4 # number of algae

    if N==4:
        ind1=np.argmin(np.array(Qps[0:100]))
        ind2=np.argmin(np.array(Qps[100:200]))
        ind3=np.argmin(np.array(Qps[200:300]))
        ind4=np.argmin(np.array(Qps[300:400]))
        #quotas
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

        # half saturation constants
        Nc1=Ncs[0:100][ind1]
        Nc2=Ncs[100:200][ind2]
        Nc3=Ncs[200:300][ind3]
        Nc4=Ncs[300:400][ind4]
        Ncs=[Nc1, Nc2, Nc3, Nc4]

        # vcell olumes
        V1=Vols[0:100][ind1]
        V2=Vols[100:200][ind2]
        V3=Vols[200:300][ind3]
        V4=Vols[300:400][ind4]
        Vs=[V1, V2, V3, V4]

    # choose phytoplankton type
    typePhytos=['Diatom', 'Eukaryote', 'Synechococcus', 'Prochlorochoccus']
    typePhyto=typePhytos[indice]

    beta=betas[indice]
    lp=lpes[indice]
    mu=mu_max[indice]
    print(mu_max)
    print(lp)
    print(beta)
    Nc=Ncs[indice]
    print(Nc)

    n_env=3

    tauT = 0.8
    bT=4
    AT=4000
    BT=0.0003
    TN=293.15

    Temps=[10, 20, 25]
    Rs=[5, 1, 0.1]
    T_deps=[tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN)) for Temp in Temps]
    R_deps=[R/(R+Nc) for R in Rs]
    
    mus=[mu*T_deps[i]*R_deps[i] for i in range(n_env)]
    
    d=0.1 # phyto mortality rate
    mui=0
    Qv=Qvs[indice]
    Qp=Qps[indice]
    print(Qp)
    print(Qv)
    eps_r=1e-6 # conversion to resistant type
    eps_lr=1e-6 # loss of resitance

    N_res=10 
    dNs=[0.5, 0.05, 0.01]
    # grazing parameters
    Sc=0.226 #dutkiewicz 2020 (1.5 in umolC.L-1)
    eps_z=0.3
    dz=0.067
    phiz=9.8

    # exploring parameters:
    # adsorption rate
    phis=list(np.logspace(-12, -8, (12-8)*9+1))
    # infection efficiency
    eps=1 # for first case
    epss=list(np.logspace(0, 4, 4*9+1))
    epss=epss[0:len(epss)-1]

    #loss of infectivity
    epso=0
    # virus linear mortality
    ms=[i/10 for i in range(1, 11)]
    ms=[0.01, 0.05]+ms
    if POM=='1':
        ms=[0.01]
    # virus quadratic mortality
    m2s_high=np.linspace(100, 2000, 20)
    m2s_low=np.linspace(10, 90, 9)
    m2s=np.concatenate((m2s_low, m2s_high))
    m2s=np.concatenate((np.zeros(1), m2s))
    # zoop quadratic mortality
    dz2s=np.linspace(2, 30, 15)
    #dz2s=np.linspace(10, 30, 21)
    dz2s[10]=22.4 # DARWIN
    # resistance strength
    res_ratios=list(np.logspace(0, 4, 4*9+1))
    res_ratios=epss[0:len(res_ratios)-1]
    res_ratios.append(0) # full resistance
    # resisiance cost
    cost_res=[i/10 for i in range(5, 10)]
    cost_res.append(0.95)
    cost_res.append(1)

    #POMc=[i/10 for i in range(1,11)]
    POMc=[i for i in range(2, 16,1)]
    POMc=POMc+[20]


    alpha=1

    # searcg grid
    if POM=='0':
        if type_m=='SIVZ':
            param_space=[(phi, m, m2, dz2, cost) for phi in phis for m in ms for m2 in m2s for dz2 in dz2s for cost in cost_res]
        elif type_m=='SIVZ_intra':
            param_space=[(phi, eps, m, m2, dz2, cost) for phi in phis for eps in epss for m in ms for m2 in m2s for dz2 in dz2s for cost in cost_res]
        elif type_m in ['SIVRZ', 'SIVRZ_intra']:    
            param_space=[(phi, m, m2, dz2, res, cost) for phi in phis for m in ms for m2 in m2s for dz2 in dz2s for res in res_ratios for cost in cost_res]
    elif POM=='1':
        if type_m=='SIVZ':
            param_space=[(phi, m, m2, dz2, cost, Pc) for phi in phis for m in ms for m2 in m2s for dz2 in dz2s for cost in cost_res for Pc in POMc]
        elif type_m=='SIVZ_intra':
            param_space=[(phi, eps, m, m2, dz2, cost, Pc) for phi in phis for eps in epss for m in ms for m2 in m2s for dz2 in dz2s for cost in cost_res for Pc in POMc]
        elif type_m in ['SIVRZ', 'SIVRZ_intra']:
            param_space=[(phi, m, m2, dz2, res, cost, Pc) for phi in phis for m in ms for m2 in m2s for dz2 in dz2s for res in res_ratios for cost in cost_res for Pc in POMc]
    
    # realistic concentration ranges and targtes
    A_cond_low_o, A_cond_high_o, V_cond_low_o, V_cond_high_o, Z_cond_low_o, Z_cond_high_o, I_cond_high_o, I_cond_low_o, perc_cond_high_o, perc_cond_low_o, target_conc_o=concentration_ranges(indice, 'oligotrophic')
    A_cond_low_m, A_cond_high_m, V_cond_low_m, V_cond_high_m, Z_cond_low_m, Z_cond_high_m, I_cond_high_m, I_cond_low_m, perc_cond_high_m, perc_cond_low_m, target_conc_m=concentration_ranges(indice, 'mesotrophic')
    A_cond_low_u, A_cond_high_u, V_cond_low_u, V_cond_high_u, Z_cond_low_u, Z_cond_high_u, I_cond_high_u, I_cond_low_u, perc_cond_high_u, perc_cond_low_u, target_conc_u=concentration_ranges(indice, 'upwelling')

    # column names of output file and output file
    if type_m=='SIVZ':
        if POM=='0':
            column_names=['valid_envs','phi', 'dv', 'dv2', 'dz2','cost', 'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u', 'ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            write_vector(column_names,'SIVZ_extracellular_res_optimization_header.txt', sep=' ')
            file_name='results_optimization_params/SIVZ_extracellular_res_optimization_'+typePhyto+'_'+chunk+'.txt'
        elif POM=='1':
            column_names=['valid_envs','phi', 'dv', 'dv2', 'dz2','cost','Pc', 'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u', 'ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            write_vector(column_names,'SIVZ_extracellular_res_optimization_header_POM.txt', sep=' ')
            file_name='results_optimization_params/SIVZ_extracellular_res_optimization_'+typePhyto+'_'+chunk+'_POM.txt'
    elif type_m=='SIVZ_intra':
        if POM=='0':
            column_names=['valid_envs','phi','eps', 'dv', 'dv2', 'dz2','cost' ,'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u','ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            write_vector(column_names,'SIVZ_intracellular_res_optimization_header.txt', sep=' ')
            file_name='results_optimization_params/SIVZ_intracellular_res_optimization_'+typePhyto+'_'+chunk+'.txt'
        elif POM=='1':
            column_names=['valid_envs','phi','eps', 'dv', 'dv2', 'dz2','cost' ,'Pc', 'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u','ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            write_vector(column_names,'SIVZ_intracellular_res_optimization_header_POM.txt', sep=' ')
            file_name='results_optimization_params/SIVZ_intracellular_res_optimization_'+typePhyto+'_'+chunk+'_POM.txt'
    elif type_m in ['SIVRZ', 'SIVRZ_intra']:
        if  POM=='0':
            column_names=['valid_envs','phi', 'dv', 'dv2', 'dz2', 'res', 'cost', 'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u','ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            if type_m=='SIVRZ':
                write_vector(column_names,'SIVRZ_extracellular_res_optimization_header.txt', sep=' ')
                file_name='results_optimization_params/SIVRZ_extracellular_res_optimization_'+typePhyto+'_'+chunk+'.txt'
            elif type_m=='SIVRZ_intra':
                write_vector(column_names,'SIVRZ_intracellular_res_optimization_header.txt', sep=' ')
                file_name='results_optimization_params/SIVRZ_intracellular_res_optimization_'+typePhyto+'_'+chunk+'.txt'
        elif POM=='1':
            column_names=['valid_envs','phi', 'dv', 'dv2', 'dz2', 'res', 'cost','Pc', 'ED_tot','ED_u', 'ED_m','ED_o','ED_av','ER_u','ER_m', 'ER_o' ,'ER_av','P_u','V_u', 'Z_u', 'Inf_u','PK_u','P_m','V_m', 'Z_m', 'Inf_m','PK_m','P_o','V_o', 'Z_o', 'Inf_o','PK_o']
            if type_m=='SIVRZ':
                write_vector(column_names,'SIVRZ_extracellular_res_optimization_header_POM.txt', sep=' ')
                file_name='results_optimization_params/SIVRZ_extracellular_res_optimization_'+typePhyto+'_'+chunk+'_POM.txt'
            elif type_m=='SIVRZ_intra':
                write_vector(column_names,'SIVRZ_intracellular_res_optimization_header_POM.txt', sep=' ')
                file_name='results_optimization_params/SIVRZ_intracellular_res_optimization_'+typePhyto+'_'+chunk+'_POM.txt'
    
    # if the output file exists, remove it
    if os.path.exists(file_name):
        os.remove(file_name)  # Remove the file
        print(f"{file_name} has been removed.")

    n_update=100 # to calculate faster the equilibrium
    follow_file='follow_optim_'+typePhyto+'_'+type_m+ '_'+POM+'_'+chunk+'.txt'
    co=0

    ncols=len(column_names)

    # parameter space corresponding to the chunk (subparameter space explored)
    n_chunks=50
    chunk_size=round(len(param_space)/n_chunks)
    if chunk==str(n_chunks):
        i1=chunk_size*(int(chunk)-1)
        i2=len(param_space)
    else:
        i1=chunk_size*(int(chunk)-1)
        i2=chunk_size*(int(chunk)-1)+chunk_size+1
    param_space=param_space[i1:i2]
    
    leno=len(phis)*len(epss)*len(m2s)*len(ms)*len(cost_res)*len(dz2s)
    print(chunk_size)

    # main lopp through the parameter space expored
    for j,k in enumerate(param_space):
        if type_m=='SIVZ':
            phi=k[0]
            m=k[1]
            m2=k[2]
            dz2=k[3]
            cost=k[4]
            if POM=='1':
                Pc=k[5]
        elif type_m=='SIVZ_intra':
            phi=k[0]
            eps=1/k[1]
            m=k[2]
            m2=k[3]
            dz2=k[4]
            cost=k[5]
            if POM=='1':
                Pc=k[6]
        elif type_m in ['SIVRZ', 'SIVRZ_intra']:
            phi=k[0]
            m=k[1]
            m2=k[2]
            dz2=k[3]
            res=k[4]
            cost=k[5]
            if POM=='1':
                Pc=k[6]
            mur=mu*cost
            if type_m=='SIVRZ':
                phir=phi/res
                epsr=1
            else:
                if res==0:
                    epsr=0
                else:
                    epsr=1/res
                phir=phi


        if POM=='0':
            phi_o=phi
            phi_m=phi
            phi_u=phi
        elif POM=='1':
            phi_o=phi*1/(1+pow(Rs[2]*2/Pc,2))
            phi_m=phi*1/(1+pow(Rs[1]*2/Pc,2))
            phi_u=phi*1/(1+pow(Rs[0]*2/Pc,2))
            if type_m=='SIVRZ':
                phir_u=phi_u/res
                phir_m=phi_m/res
                phir_o=phi_o/res
        
        # variables dependent on growth rate
        KCs=[(-dNs[i]*Rs[i]+dNs[i]*N_res)/(mus[i]*cost) for i in range(n_env)]
        CCs=[KCs[i]/(mus[i]*cost-d*T_deps[i]) for i in range(n_env)]
        S_deps=[pow(KCs[i],2)/(pow(KCs[i],2)+pow(Sc, 2)) for i in range(n_env)]
        phizs=[phiz*T_deps[i]*S_deps[i] for i in range(n_env)]
        
        # equilibrium for the different environments
        if type_m in ['SIVZ', 'SIVZ_intra']:
            if m2!=0:
                S_star_u, I_star_u, V_star_u, Z_star_u, case_u=equilibrium_SIVZ_m2(mus[0]*cost, mui, lp, beta, phi_u, d*T_deps[0], m*T_deps[0],m2*T_deps[0], Qv, Qp,  eps, epso,phizs[0],eps_z,dz*T_deps[0],dz2*T_deps[0], CCs[0], n_update)
                S_star_m, I_star_m, V_star_m, Z_star_m, case_m=equilibrium_SIVZ_m2(mus[1]*cost, mui, lp, beta, phi_m, d*T_deps[1], m*T_deps[1],m2*T_deps[1], Qv, Qp,  eps, epso,phizs[1],eps_z,dz*T_deps[1],dz2*T_deps[1], CCs[1], n_update)
                S_star_o, I_star_o, V_star_o, Z_star_o, case_o=equilibrium_SIVZ_m2(mus[2]*cost, mui, lp, beta, phi_o, d*T_deps[2], m*T_deps[2],m2*T_deps[2], Qv, Qp,  eps, epso,phizs[2],eps_z,dz*T_deps[2],dz2*T_deps[2], CCs[2], n_update)
            else:
                S_star_u, I_star_u, V_star_u, Z_star_u=equilibrium_SIVZ(mus[0]*cost, mui, lp, beta, phi_u, d*T_deps[0], m*T_deps[0],m2*T_deps[0], Qv, Qp,Qz,  eps, epso,phizs[0],eps_z,dz*T_deps[0],dz2*T_deps[0], CCs[0])
                if phi_u*KCs[0]/Qp<10*d*T_deps[0]: 
                    case_u=1
                else:
                    case_u=0
                S_star_m, I_star_m, V_star_m, Z_star_m=equilibrium_SIVZ(mus[1]*cost, mui, lp, beta, phi_m, d*T_deps[1], m*T_deps[1],m2*T_deps[1], Qv, Qp,Qz,  eps, epso,phizs[1],eps_z,dz*T_deps[1],dz2*T_deps[1], CCs[1])
                if phi_m*KCs[1]/Qp<10*d*T_deps[1]:
                    case_m=1
                else:
                    case_m=0
                S_star_o, I_star_o, V_star_o, Z_star_o=equilibrium_SIVZ(mus[2]*cost, mui, lp, beta, phi_o, d*T_deps[2], m*T_deps[2],m2*T_deps[2], Qv, Qp,Qz,  eps, epso,phizs[2],eps_z,dz*T_deps[2],dz2*T_deps[2], CCs[2])
                if phi_o*KCs[2]/Qp<10*d*T_deps[2]:
                    case_o=1
                else:
                    case_o=0

            R_star_o, R_star_m, R_star_u=0,0,0
            c_u=S_star_u>0 and I_star_u>0 and V_star_u>0 and Z_star_u> 0 and case_u==1
            c_m=S_star_m>0 and I_star_m>0 and V_star_m>0 and Z_star_m> 0 and case_m==1
            c_o= S_star_o>0 and I_star_o>0 and V_star_o>0 and Z_star_o>0 and case_o==1
            c_val=(c_m and c_o)
        elif type_m in ['SIVRZ', 'SIVRZ_intra']:
            if m2!=0:
                S_star_u, I_star_u, V_star_u, R_star_u, Z_star_u, case_u=equilibrium1_SIVRZ_m2(mus[0], mui, lp, beta, phi_u, d*T_deps[0], m*T_deps[0],m2*T_deps[0], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_u,mur,phizs[0],eps_z,dz*T_deps[0],dz2*T_deps[0], CCs[0], alpha)
                S_star_m, I_star_m, V_star_m, R_star_m, Z_star_m, case_m=equilibrium1_SIVRZ_m2(mus[1], mui, lp, beta, phi_m, d*T_deps[1], m*T_deps[1],m2*T_deps[1], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_m,mur,phizs[1],eps_z,dz*T_deps[1],dz2*T_deps[1], CCs[1], alpha)
                S_star_o, I_star_o, V_star_o, R_star_o, Z_star_o, case_o=equilibrium1_SIVRZ_m2(mus[2], mui, lp, beta, phi_o, d*T_deps[2], m*T_deps[2],m2*T_deps[2], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_o,mur,phizs[2],eps_z,dz*T_deps[2],dz2*T_deps[2], CCs[2], alpha)
            else:
                S_star_u, I_star_u, V_star_u, R_star_u, Z_star_u=equilibrium1_SIVRZ(mus[0], mui, lp, beta, phi_u, d*T_deps[0], m*T_deps[0],m2*T_deps[0], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_u,mur,phizs[0],eps_z,dz*T_deps[0],dz2*T_deps[0], CCs[0], alpha)
                if phi_u*KCs[0]/Qp<10*d*T_deps[0]:
                    case_u=1
                else:
                    case_u=0
                S_star_m, I_star_m, V_star_m, R_star_m, Z_star_m=equilibrium1_SIVRZ(mus[1], mui, lp, beta, phi_m, d*T_deps[1], m*T_deps[1],m2*T_deps[1], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_m,mur,phizs[1],eps_z,dz*T_deps[1],dz2*T_deps[1], CCs[1], alpha)
                if phi_m*KCs[1]/Qp<10*d*T_deps[1]:
                    case_m=1
                else:
                    case_m=0
                S_star_o, I_star_o, V_star_o, R_star_o, Z_star_o=equilibrium1_SIVRZ(mus[2], mui, lp, beta, phi_o, d*T_deps[2], m*T_deps[2],m2*T_deps[2], Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir_o,mur,phizs[2],eps_z,dz*T_deps[2],dz2*T_deps[2], CCs[2], alpha)
                if phi_o*KCs[2]/Qp<10*d*T_deps[2]:
                    case_o=1
                else:
                    case_o=0
            c_u=S_star_u>0 and I_star_u>0 and V_star_u>0 and Z_star_u> 0 and R_star_u>0 and case_u==1
            c_m=S_star_m>0 and I_star_m>0 and V_star_m>0 and Z_star_m> 0 and R_star_m>0 and case_m==1
            c_o= S_star_o>0 and I_star_o>0 and V_star_o>0 and Z_star_o>0 and R_star_o>0 and case_o==1
            c_val=(c_m and c_o) 

        # check positivity of equilibrium in both the mesotrophic and oligotrophic environment
        if c_val:
            v1=target_conc_o+target_conc_m+target_conc_u
            v2=[(S_star_o+I_star_o+R_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz, (S_star_m+I_star_m+R_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz, (S_star_u+I_star_u+R_star_u)/Qp, V_star_u/Qv, Z_star_u/Qz]
            
            #euclidean distance for each env and tot distance
            ed_tot=euclidean_distance(v1, v2)
            
            ed_u=euclidean_distance(v1[6:9], v2[6:9])
            ed_m=euclidean_distance(v1[3:6], v2[3:6])
            ed_o=euclidean_distance(v1[0:3], v2[0:3])

            if indice in [2,3]:
                ed_av=(ed_m+ed_o)/2
            elif indice in [0,1]:
                ed_av=(2*ed_m+ed_o)/3
                
            # absolute error to target for each env 
            er_u=absolute_error(v1[6:9], v2[6:9])
            er_m=absolute_error(v1[3:6], v2[3:6])
            er_o=absolute_error(v1[0:3], v2[0:3])

            # consider only oligotrophic and mesotrophic env for average distance, give more weitgh to mesotrophic environment for diatom and eukaryote
            if indice in [2,3]:
                er_av=(er_m+er_o)/2
            elif indice in [0,1]:
                er_av=(2*er_m+er_o)/3

            perc_inf_o=I_star_o*100/(I_star_o+S_star_o+R_star_o)
            perc_inf_m=I_star_m*100/(I_star_m+S_star_m+R_star_m)
            perc_inf_u=I_star_u*100/(I_star_u+S_star_u+R_star_u)

            mZ_o=phizs[2]*Z_star_o
            mZ_m=phizs[1]*Z_star_m
            mZ_u=phizs[0]*Z_star_u
        
            mV_o=perc_inf_o/100/lp
            mV_m=perc_inf_m/100/lp
            mV_u=perc_inf_u/100/lp

            perc_kill_o=mV_o/(mV_o+mZ_o+d*T_deps[2])
            perc_kill_m=mV_m/(mV_m+mZ_m+d*T_deps[1])
            perc_kill_u=mV_u/(mV_u+mZ_u+d*T_deps[0])

        
            # chack for each env if equilibrium concentration are in the range of realistic concentration 
            val_u=0
            if (S_star_u+I_star_u+R_star_u)/Qp > A_cond_low_u and (S_star_u+I_star_u+R_star_u)/Qp<A_cond_high_u and V_star_u/Qv>V_cond_low_u and V_star_u/Qv<V_cond_high_u and Z_star_u/Qz>Z_cond_low_u and Z_star_u/Qz<Z_cond_high_u and perc_inf_u>I_cond_low_u and perc_inf_u<I_cond_high_u and perc_kill_u>perc_cond_low_u and perc_kill_u<perc_cond_high_u:
                val_u=1
            val_m=0
            if  (S_star_m+I_star_m+R_star_m)/Qp > A_cond_low_m and (S_star_m+I_star_m+R_star_m)/Qp<A_cond_high_m and V_star_m/Qv>V_cond_low_m and V_star_m/Qv<V_cond_high_m and Z_star_m/Qz>Z_cond_low_m and Z_star_m/Qz<Z_cond_high_m and perc_inf_m>I_cond_low_m and perc_inf_m<I_cond_high_m and perc_kill_m>perc_cond_low_m and perc_kill_m<perc_cond_high_m:
                val_m=1
            val_o=0
            if  (S_star_o+I_star_o+R_star_o)/Qp > A_cond_low_o and (S_star_o+I_star_o+R_star_o)/Qp<A_cond_high_o and V_star_o/Qv>V_cond_low_o and V_star_o/Qv<V_cond_high_o and Z_star_o/Qz>Z_cond_low_o and Z_star_o/Qz<Z_cond_high_o and perc_inf_o>I_cond_low_o and perc_inf_o<I_cond_high_o and perc_kill_o>perc_cond_low_o and perc_kill_o<perc_cond_high_o:
                val_o=1

            val_t=val_m+val_o
            valid_envs=str(val_u)+str(val_m)+str(val_o)
            valid_eq=0

            # valid if oligotrophic and mesotrophic equilibria are both in the realistic range
            if val_t==2:
                valid_eq=1

            # write to file
            if valid_eq==1:
                if type_m=='SIVZ':
                    if POM=='0':
                        data_to_write=[valid_envs, phi, m, m2, dz2, cost,ed_tot, ed_u, ed_m, ed_o, ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                    elif POM=='1':
                        data_to_write=[valid_envs, phi, m, m2, dz2, cost, Pc,ed_tot, ed_u, ed_m, ed_o, ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                elif type_m=='SIVZ_intra':
                    if POM=='0':
                        data_to_write=[valid_envs,phi, eps, m, m2, dz2, cost,ed_tot, ed_u, ed_m, ed_o, ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                    elif POM=='1':
                        data_to_write=[valid_envs,phi, eps, m, m2, dz2, cost, Pc,ed_tot, ed_u, ed_m, ed_o, ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                elif type_m in ['SIVRZ', 'SIVRZ_intra']:
                    if POM=='0':
                        data_to_write=[valid_envs, phi, m, m2, dz2, res,cost,ed_tot, ed_u, ed_m, ed_o,ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                    elif POM=='1':
                        data_to_write=[valid_envs, phi, m, m2, dz2, res,cost, Pc,ed_tot, ed_u, ed_m, ed_o,ed_av, er_u, er_m, er_o, er_av,(S_star_u+I_star_u)/Qp,  V_star_u/Qv, Z_star_u/Qz,perc_inf_u,perc_kill_u , (S_star_m+I_star_m)/Qp, V_star_m/Qv, Z_star_m/Qz,perc_inf_m,perc_kill_m, (S_star_o+I_star_o)/Qp, V_star_o/Qv, Z_star_o/Qz,perc_inf_o,perc_kill_o]
                if len(data_to_write)==ncols:
                    with open(file_name, 'a') as f:
                        if co==0:
                            for fl in data_to_write:
                                f.write(str(fl))
                                f.write(' ')
                                co=co+1
                            f.write('\n')
                        else:
                            for fl in data_to_write:
                                f.write(str(fl))
                                f.write(' ')
                            f.write('\n')
