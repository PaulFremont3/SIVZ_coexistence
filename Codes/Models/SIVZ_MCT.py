import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import copy
from numpy import inf
import random as random
from SIVZ_functions import *


if __name__ == '__main__':
    indice=int(sys.argv[1]) # phytoplankton type, 0: small diatom, 1, picoeukaryote, 2, synechococcus, 3, prochlorococcus
    eff=sys.argv[2] # loss of infection rate (0)
    otype=sys.argv[3] # ocean type
    dz2=float(sys.argv[4]) # quadratic mortality term of zooplankton
    alph=float(sys.argv[5])  # strength of intraguild predation 0-1
    dv2=float(sys.argv[6])  # quadratic mortality term of virus

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
    # zooplankton nitrogen quota
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

    # choose the 4 smallest cells of each type (diatom, euk, syn and pro)
    N=4 # number of algae
    if N==4:
        ind1=np.argmin(np.array(Qps[0:100]))
        ind2=np.argmin(np.array(Qps[100:200]))
        ind3=np.argmin(np.array(Qps[200:300]))
        ind4=np.argmin(np.array(Qps[300:400]))
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

        # maximum growth rate
        mu1=mu_max[0:100][ind1]
        mu2=mu_max[100:200][ind2]
        mu3=mu_max[200:300][ind3]
        mu4=mu_max[300:400][ind4]
        mu_max=[mu1,mu2,mu3,mu4]

        # half saturation constant for nutrient
        Nc1=Ncs[0:100][ind1]
        Nc2=Ncs[100:200][ind2]
        Nc3=Ncs[200:300][ind3]
        Nc4=Ncs[300:400][ind4]
        Ncs=[Nc1, Nc2, Nc3, Nc4]

    print(mu_max)

    # selct phytpplankton type
    typePhytos=['Diatom', 'Eukaryote', 'Synechococcus', 'Prochlorochoccus']
    typePhyto=typePhytos[indice]

    bs=betas[indice]

    
    suffix=typePhyto+'_BS'+str(bs)+'_LOI'+eff
    if dz2==0:
        suffix+='_no-dz2'
    if dv2==0:
        suffix+='_no-dv2'
    else:
        suffix+='dv2-'+str(dv2)
    if alph!=1:
        suffix+='_IG'+str(alph)

    # Temperature dependency
    T_dep=1 # no T_dep in idealized case
    R=1 # default nutrient
    if otype != '0':
        suffix+='_'+otype
        tauT = 0.8
        bT=4
        AT=4000
        BT=0.0003
        TN=293.15
        if otype=='upwelling':
            Temp=10
            R=5
            T_dep=tauT*np.exp(-AT*(1/(273.15+Temp)-1/TN))
            mu_max=[mu_max[i]*T_dep*R/(R+Ncs[i]) for i in range(N)]
            print(mu_max)
        elif otype=='oligotrophic':
            Temp=30
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
    # simulation params
    dt=1/(48)
    nyears=20
    ndays=round(365*nyears)

    mu=mu_max[indice]
    mui=0 # infected phyto max growth rate
    beta=float(bs) # burst size
    d=0.1*T_dep # phyto mortality rate
    m=0.1*T_dep #virus mortality
    m2=dv2*T_dep #virus qudratic mortality
    Qv=Qvs[indice]
    Qp=Qps[indice]
    Nc=Ncs[indice]
    eps=1 # adsorption efficiency coefficient
    epso=float(eff) # loss of infectivity rate
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
    KC_s=(-dN*R+SN)*(R+Nc)/(mu*R)
    CC=KC_s/(mu-d)
 
    print(KC_s)
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

    # adsorption rate parameter space
    phis=list(np.logspace(-12, -8, (12-8)*9+1))

    #latent period parameter space
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]

    # paraemter space grid
    grid_test = [(ph, lp) for ph in phis for lp in lps]

    # matrice to store effetcs from MCTs
    ntot=6
    matrices_effects_v_inv=[np.zeros( (len(phis),len(lps)) ) for i in range(ntot)]
    matrices_effects_z_inv=[np.zeros( (len(phis),len(lps)) ) for i in range(ntot)]
    matrices_effects_stab=[np.zeros( (len(phis),len(lps)) ) for i in range(ntot)]
    invasibilty_both=np.zeros( (len(phis),len(lps)) )
    invasibilty_both_binary=np.zeros( (len(phis),len(lps)) )

    # main loop trhough parameter space
    for j, k in enumerate(grid_test):
        phi=k[0]
        lp=k[1]
        
        i=np.where(np.array(phis)==phi)[0][0]
        j=np.where(np.array(lps)==k[1])[0][0]

        # Modern Coexistence Theory analysis
        result=MCT_analysis_SIVZ(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, ntot)

        # fill matrices
        for l in range(ntot):
            if l==0:
                matrices_effects_v_inv[l][i,j]=sum(result[0:(ntot-1)])
                matrices_effects_z_inv[l][i,j]=sum(result[(ntot-1):(2*ntot-3)])
                matrices_effects_stab[l][i,j]=( sum(result[0:(ntot-1)]) + sum(result[(ntot-1):(2*ntot-3)]) )/2
                if sum(result[0:(ntot-1)]) >0  and sum(result[(ntot-1):(2*ntot-3)]) >0:
                    invasibilty_both[i,j] = sum(result[0:(ntot-1)]) + sum(result[(ntot-1):(2*ntot-3)])
                    invasibilty_both_binary[i,j]=1
                else:
                    invasibilty_both[i,j]=np.nan
                    if np.isnan(sum(result[0:(ntot-1)])) or np.isnan(sum(result[(ntot-1):(2*ntot-3)])):
                        invasibilty_both_binary[i,j]=np.nan
            if l>0:
                matrices_effects_v_inv[l][i,j]=result[l-1]
                matrices_effects_z_inv[l][i,j]=result[ntot+l-2]
                matrices_effects_stab[l][i,j]=np.mean(np.array([result[l-1], result[ntot+l-2]]))
        
    # pdf eith effects over the paraneter space
    pp = PdfPages('SIVZ_model_phi_latent_period_'+suffix+'_MCT.pdf')
    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]
    

    mxs=[]
    mis=[]
    for l in range(ntot):
        mat1=matrices_effects_z_inv[l]
        mat2=matrices_effects_v_inv[l]
        mat1[mat1<-100]=-100
        mat2[mat2<-100]=-100
        mat1[mat1>100]=100
        mat2[mat2>100]=100
        mxs.append(np.max(mat1))
        mxs.append(np.max(mat2))

        mat3=matrices_effects_stab[l]
        mat3[mat3<-100]=-100
        mat3[mat3>100]=100
        mxs.append(np.max(mat3))

        mis.append(np.min(mat1))
        mis.append(np.min(mat2))
        mis.append(np.min(mat3))

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

    titles=['Invasion growth rate','Fluctuation free growth rate', 'Relative non linearity in S', 'Relative non linearity in I', 'I, S covariance', 'I, S variance interaction']
    for l in range(ntot):
        mat=matrices_effects_z_inv[l]
        mat[mat<-100]=-100
        mat[mat>100]=100
        maxi=np.nanmax(mat)
        mini=np.nanmin(mat)
        amx=max([abs(mini), abs(maxi)])
        ma=amx
        nma=-amx
        if (maxi==mini):
            ma=0.01
            nma=-0.01
        bounds=np.linspace(-amx, amx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        plot_with_scale(mat,cmap0,nma,ma, atickx, aticky, alabel_tickx, alabel_ticky, titles[l]+' (Z invading)', norm=norm)
        pp.savefig()
    for l in range(ntot):
        mat=matrices_effects_v_inv[l]
        mat[mat<-100]=-100
        mat[mat>100]=100
        maxi=np.nanmax(mat)
        mini=np.nanmin(mat)
        amx=max([abs(mini), abs(maxi)])
        ma=amx
        nma=-amx
        if (maxi==mini):
            ma=0.01
            nma=-0.01
        bounds=np.linspace(-amx, amx, 100)
        print(titles[l])
        print('max')
        print(maxi)
        print(mini)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        plot_with_scale(mat,cmap0,nma,ma, atickx, aticky, alabel_tickx, alabel_ticky, titles[l]+' (V invading)', norm=norm)
        pp.savefig()
    for l in range(ntot):
        mat=matrices_effects_stab[l]
        mat[mat<-100]=-100
        mat[mat>100]=100
        maxi=np.nanmax(mat)
        mini=np.nanmin(mat)
        amx=max([abs(mini), abs(maxi)])
        ma=amx
        nma=-amx
        if (maxi==mini):
            ma=0.01
            nma=-0.01
        bounds=np.linspace(-amx, amx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        plot_with_scale(mat,cmap0,nma,ma, atickx, aticky, alabel_tickx, alabel_ticky, titles[l]+' (stab)', norm=norm)
        pp.savefig()

    cmap = matplotlib.colormaps['PuOr_r']
    colors0 = cmap(np.linspace(0, 1, 200))
    colors0=colors0[100:]
    cmap0=matplotlib.colors.ListedColormap(colors0)

    amx=np.nanmax(np.array(invasibilty_both))
    bounds=np.linspace(0, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=0
    ma=amx
    plot_with_scale(invasibilty_both, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Total invasion growth rate (>0)', norm=norm)
    pp.savefig()

    mi=0
    mx=1
    cmap = matplotlib.colors.ListedColormap(['white', 'green'])
    plot_with_scale(invasibilty_both_binary, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Invasibility', norm=norm)
    pp.savefig()

    pp.close()
    
        


