import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import copy
from numpy import inf
from SVZ_functions import *
from generic_functions import *

#if __name__ == '__main__':
def main():
    print("=== Script started ===")

    indice=int(sys.argv[1]) # choose phytoplantkon: 0, small diatom, 1, picoeukaryote, 2 synechococcus, 3, prochlorococcus
    eff=sys.argv[2] # loss of infectivity (irrelevant here)
    otype=sys.argv[3] # ocean type 
    dz2=float(sys.argv[4]) # quadratic mortality of zooplankton
    m2=float(sys.argv[5]) # quadratic mortality of virus

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
    # zooplankton quota
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
    # choose life history traits of the smallest of each type of phytoplankton
    if N==4:
        ind1=np.argmin(np.array(Qps[0:100]))
        ind2=np.argmin(np.array(Qps[100:200]))
        ind3=np.argmin(np.array(Qps[200:300]))
        ind4=np.argmin(np.array(Qps[300:400]))
        # nitrogen quots
        Qp1=Qps[0:100][ind1]
        Qp2=Qps[100:200][ind2]
        Qp3=Qps[200:300][ind3]
        Qp4=Qps[300:400][ind4]
        # burst sizes
        bs1=betas[0:100][ind1]
        bs2=betas[100:200][ind2]
        bs3=betas[200:300][ind3]
        bs4=betas[300:400][ind4]
        betas=[bs1,bs2,bs3,bs4]
        # latent periods
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
        # volumes
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
    R=1 # default nutrient concentration
    T_dep=1 # no T_dep in idealized case
    # suffix to identify simulation and output pdf
    suffix=typePhyto+'_BS'+str(bs)+'_LOI'+eff
    if dz2==0:
        suffix+='_no-dz2'
    if m2!=0:
        suffix+='_m2-'+str(m2)
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
    nyears=25
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

    # grazing parameters
    phiz=9.8*T_dep
    eps_z=0.3
    dz=0.067*T_dep
    dz2=dz2*T_dep

    print('Carying capacity')
    print(KC_s)

    if indice in [1,2,3]:
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    elif indice in [0]:
        phis=list(np.logspace(-11, -7, (12-8)*9+1))

    # no latent period
    lps=[0]

    # parameter space tested
    grid_test = [(ph, lp) for ph in phis for lp in lps]

    alphas=[0, 0.1, 0.5, 1]

    n_state=5
    # matrices to store results
    limit_val_mat = np.zeros( (len(phis),len(lps)) )
    mortality_ratios = np.zeros( (len(phis),len(lps)) )

    a_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    o_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    f_stable_states_matrix = [np.empty( (len(phis),len(lps)) ) for i in range(n_state)]
    for l in range(n_state):
        a_stable_states_matrix[l][:]=np.nan
        o_stable_states_matrix[l][:]=np.nan
        f_stable_states_matrix[l][:]=np.nan

    final_S = np.zeros( (len(phis),len(lps)) )
    final_Surv = np.zeros( (len(phis),len(lps)) )
    final_Surv_alpha = [np.zeros( (len(phis),len(lps)) ) for i in range(len(alphas))]
    final_ZV_ratio = np.zeros( (len(phis),len(lps)) )
    final_ZV_ratio_count = np.zeros( (len(phis),len(lps)) )
    final_Z=np.zeros( (len(phis),len(lps)) )
    final_V=np.zeros( (len(phis),len(lps)) )

    exclusion_times=np.zeros( (len(phis),len(lps)) )
    final_pers=np.zeros( (len(phis),len(lps)) )
    final_mods=np.zeros( (len(phis),len(lps)) )
    final_stab=np.empty( (len(phis),len(lps)) )
    final_stab[:]=np.nan
    final_unstab_lps=[]
    final_unstab_phis=[]
    theor_surv_phi=[]
    theor_surv_lp=[]
    theor_surv_phi_10=[]
    theor_surv_lp_10=[]

    # concentration ranges
    A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc=concentration_ranges(indice, otype)

    valid_concentrations=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_bis=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_ter=np.zeros( (len(phis),len(lps)) )
    valid_concentrations_4=np.zeros( (len(phis),len(lps)) )

    distance_to_target=np.zeros( (len(phis),len(lps)) )
    distance_to_target_A=np.zeros( (len(phis),len(lps)) )
    distance_to_target_V=np.zeros( (len(phis),len(lps)) )
    distance_to_target_Z=np.zeros( (len(phis),len(lps)) )

    # main loop: explore parameter space
    init_conditions=[0.001,0.001, 0.001]
    for k in grid_test:
        phi=k[0]
        lp=k[1]

        i=np.where(np.array(phis)==phi)[0][0]
        j=np.where(np.array(lps)==k[1])[0][0]

        # simulation
        eta=lp
        result=simulation_SVZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC, dt, ndays, init_conditions)

        S=result[0]
        V=result[1]
        Z=result[2]
        mZs=result[3]
        mVs=result[4]

        # extract last year features
        i1=max([0, int(365*(nyears-1)/dt)+1])
        i2=len(V)
        V = list(map(lambda x: 0 if x == float('inf') else x, V))
        S = list(map(lambda x: 0 if x == float('inf') else x, S))
        Z = list(map(lambda x: 0 if x == float('inf') else x, Z))
        Sa=np.array(S[i1:i2])
        Va=np.array(V[i1:i2])
        Za=np.array(Z[i1:i2])

        mZa=np.array(mZs[i1:i2])
        mVa=np.array(mVs[i1:i2])

        # criterion of existence for each tracer
        surv_S=np.max(Sa/Qp)>1 
        surv_V = np.max(Va/Qv)>1 
        surv_Z=np.max(Za/Qz)>1 
        # fill matrices
        if surv_S or surv_V:
            final_S[i,j]=np.log10(np.nanmean(Sa))
        else:
            final_S[i,j]=float("NaN")
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

        # final state of the system
        if not surv_V  and not surv_Z and not surv_S:
            surv=0
        elif not surv_V  and not surv_Z  and surv_S:
            surv=1
        elif surv_V and not surv_Z:
            surv=2
        elif not surv_V and surv_Z:
            surv=3
        elif surv_V  and surv_Z:
            surv=4

        # exclusion time
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
        else:
            final_ZV_ratio_count[i,j]=float("NaN")
            mortality_ratios[i,j]=float("NaN")

        if surv==4 or surv==3:
            final_Z[i,j] = np.log10(np.nanmean(Za))
        else:
            final_Z[i,j]= float("NaN")
        if surv_V:
            mVa=np.mean(Va)
            final_V[i,j] = np.log10(mVa)
        else:
            final_V[i,j] = float("NaN")
        final_Surv[i,j]=surv

        # comparison to realistic concentration and target concentrations
        mean_A=np.nanmean(Sa)/Qp
        mean_V=np.nanmean(Va)/Qv
        mean_Z=np.nanmean(Za)/Qz
        vals_r=[mean_A, mean_V, mean_Z]
        if surv==4:

            if mean_A > A_cond_low and mean_A < A_cond_high and mean_V>V_cond_low and mean_V<V_cond_high and mean_Z>Z_cond_low and mean_Z<Z_cond_high:
                valid_concentrations[i,j]=1
            if valid_concentrations[i,j]==1 and mortality_ratios[i,j]>perc_cond_low and mortality_ratios[i,j]<perc_cond_high:
                valid_concentrations_bis[i,j]=1
            if valid_concentrations[i,j]==1 and mean_V/mean_A>=1:
                valid_concentrations_ter[i,j]=1
            if valid_concentrations[i,j]==1 and valid_concentrations_bis[i,j]==1 and valid_concentrations_ter[i,j]==1:
                valid_concentrations_4[i,j]=1

            distance_to_t = absolute_error(target_conc, vals_r)

            distance_to_target[i,j]=distance_to_t
        else:
            distance_to_target[i,j]=float("NaN")
        if surv_S or surv_V:
            distance_to_target_A[i,j]=vals_r[0]-target_conc[0]
        else:
            distance_to_target_A[i,j]=float("NaN")
        if surv_V:
            distance_to_target_V[i,j]=vals_r[1]-target_conc[1]
        else:
            distance_to_target_V[i,j]=float("NaN")
        if surv_Z:
            distance_to_target_Z[i,j]=vals_r[2]-target_conc[2]
        else:
            distance_to_target_Z[i,j]=float("NaN")

    # time series pdf
    pp = PdfPages('SVZ_model_phi_latent_period_time_series_'+suffix+'.pdf')
    if indice in [1,2,3]:
        phi_to_plot=[1e-12, 1e-11, 1e-10, 1e-9, 1e-8]
    else:
        phi_to_plot=[1e-11, 1e-10, 1e-9, 1e-8, 1e-7]
    for phi in phi_to_plot:
        tit='phi='+str(phi)
        d=0.1
        lat_per=0
        nyears0=5
        ndays=365*nyears0
        result=simulation_SVZ_rk4(mu, mui, lat_per, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,dt, ndays , init_conditions )
        S=result[0]
        V=result[1]
        Z=result[2]
        i1=0
        i2=len(S)-1
        make_plots(S, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz)
        i3=round((365)/dt)
        make_plots(S, V, Z, i1, i3,tit,dt,pp, Qp, Qv, Qz)
    pp.close()

    atickx=[i*9 if i>0 else 0 for i in range(5)]
    aticky=[]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= []

    final_S[final_S==-inf]=float("NaN")
    final_Z[final_Z==-inf]=float("NaN")
    final_V[final_V==-inf]=float("NaN")

    # main pdf with matrices plots in the explored parameter sapce
    pp = PdfPages('SVZ_model_phi_latent_period_'+suffix+'.pdf')
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
    colors[0]=(0,0,0)
    #colors[4] = (0.502,0,0)
    colors[10]=(0.8039, 0.5216, 0.2471)
    idis=[0,2,6,4,9]
    colors=[colors[idi] for idi in idis]
    colors=tuple(colors)
    cmap_R= matplotlib.colors.ListedColormap(colors)
    axi=plot_with_scale(final_Surv,cmap_R,mi,mx, atickx, aticky, alabel_tickx, alabel_ticky, 'State reached', norm=norm, yl=ylb)
    coex_mat=np.zeros( (len(phis),len(lps)) )
    coex_mat[final_Surv==4]=1
    coex_mat=np.transpose(coex_mat)
    #print(coex_mat)
    draw_borders(coex_mat, 1.5, 0.02, 'white', axi)
    pp.savefig()

    mi=np.nanmin(np.array(exclusion_times))
    mx=np.nanmax(np.array(exclusion_times))
    plot_with_scale(exclusion_times, 'YlGnBu', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Exclusion time (year)', yl=ylb)
    pp.savefig()

    exclusion_times_bis=copy.deepcopy(exclusion_times)
    exclusion_times_bis[exclusion_times_bis>1.5]=1.5
    plot_with_scale(exclusion_times_bis, 'YlGnBu', mi, 1.5, atickx, aticky, alabel_tickx, alabel_ticky , 'Exclusion time (year)', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z))
    mx=np.nanmax(np.array(final_Z))
    plot_with_scale(final_Z, 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_Z)-mt.log10(Qz))
    mx=np.nanmax(np.array(final_Z)-mt.log10(Qz))
    plot_with_scale(np.array(final_Z)-mt.log10(Qz), 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Zooplankton (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V))
    mx=np.nanmax(np.array(final_V))
    plot_with_scale(final_V, 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_V)-mt.log10(Qv))
    mx=np.nanmax(np.array(final_V)-mt.log10(Qv))
    plot_with_scale(np.array(final_V)-mt.log10(Qv), 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Virus (log10(ind.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_S))
    mx=np.nanmax(np.array(final_S))
    plot_with_scale(final_S, 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible (log10(umol.L-1))', yl=ylb)
    pp.savefig()

    mi=np.nanmin(np.array(final_S)-mt.log10(Qp))
    mx=np.nanmax(np.array(final_S)-mt.log10(Qp))
    plot_with_scale(np.array(final_S)-mt.log10(Qp), 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Susceptible (log10(ind.L-1))', yl=ylb)
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
    mx=100
    plot_with_scale_bis(mortality_ratios, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% virus induced death', yl=ylb)
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

    to_plot=[valid_concentrations, valid_concentrations_bis, valid_concentrations_ter, valid_concentrations_4]
    titls=['Valid concentrations', 'Valid concentrations + % virus kill', 'Valid concentrations + % V/A ratio', 'Valid concentrations + % virus kill+ % V/A ratio']
    mi=0
    mx=1
    cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue'])
    for k,mato in enumerate(to_plot):
        plot_with_scale(mato, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, titls[k], yl=ylb)
        pp.savefig()

    try:
        ind_min=np.nanargmin(distance_to_target)
    except:
        ind_min=np.nan
        min_d=np.nan
        phi_min=np.nan
    if not np.isnan(ind_min):
        min_d=distance_to_target[ind_min]
        min_d=min_d[0]
        phi_min=phis[ind_min]
   
        try:
            distance_to_target_val=copy.deepcopy(distance_to_target)
            distance_to_target_val[valid_concentrations_4==0]=np.nan
            ind_min_val=np.nanargmin(distance_to_target_val)
        except:
            ind_min_val,min_d_val, phi_min_val=np.nan,np.nan,np.nan
        
        if not np.isnan(ind_min_val):
            phi_min_val=phis[ind_min_val]
            min_d_val=distance_to_target[ind_min_val]
            min_d_val=min_d_val[0]
            valid_concentrations_4[ind_min_val]=2

            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations + % virus kill+ % V/A ratio', yl=ylb)
            pp.savefig()

            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_index_par]=3
            
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +% virus kill + V/A ratio', yl=ylb)
            pp.savefig()
            
            cond=(np.nan_to_num(distance_to_target_val, nan=np.inf) < min_d_val + 10) & (np.nan_to_num(valid_concentrations_4, nan=np.inf) != 3)
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
            min_d_val=np.nan
            phi_min_val=np.nan
            phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
            phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
            valid_concentrations_4[phi_index_par]=3
            valid_concentrations_4[phi_a_index_par]=4
            mi=0
            mx=np.nanmax(valid_concentrations_4)
            cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
            plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +% virus kill + V/A ratio', yl=ylb)
            pp.savefig()
    else:
        min_d_val=np.nan
        phi_min_val=np.nan
        phi_index_par=np.nanargmin(np.absolute(mt.log10(phi_e)-np.log10(np.array(phis))))
        phi_a_index_par=np.nanargmin(np.absolute(mt.log10(phi_e_a)-np.log10(np.array(phis))))
        valid_concentrations_4[phi_index_par]=3
        valid_concentrations_4[phi_a_index_par]=4
        mi=0
        mx=np.nanmax(valid_concentrations_4)
        cmap = matplotlib.colors.ListedColormap(['crimson', 'royalblue', 'gold', 'limegreen', 'cyan'])
        plot_with_scale(valid_concentrations_4, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Valid concentrations +% virus kill + V/A ratio', yl=ylb)
        pp.savefig()

    pp.close()

    if 'min_d_val' not in globals():
        min_d_val=np.nan

    print([min_d, phi_min, min_d_val, phi_min_val])
    if np.isnan(min_d):
        to_write=[min_d, phi_min, min_d_val, phi_min_val, np.nan]
    else:
        if 'vl' not in globals():
            vl=np.nan
        else:
            vl=vl[0]
        to_write=[min_d, phi_min, min_d_val, phi_min_val, vl]
    if nyears==25:
        write_vector(to_write, 'optimums_SVZ_'+suffix+'.txt', ' ')
    
    exclusion_times=np.transpose(exclusion_times)
    if nyears==25:
        write_matrix(exclusion_times, 'model_data/exclusion_times_SVZ_'+suffix+'.txt', ' ')
        write_matrix(final_Surv, 'model_data/state_reached_SVZ_'+suffix+'.txt', ' ')
    
    print("=== Script finished ===")

if __name__ == '__main__':
    main() 


