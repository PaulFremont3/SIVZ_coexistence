import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
import copy
from numpy import inf
import math as mt
from SIVZ_functions import *
from scipy.signal import find_peaks
from generic_functions import *

if __name__ == '__main__':
    mods=str(sys.argv[1]) #  main_models 

    #if mods=='main_models':
    suffixes=['SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt']
    titles=['SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic']
    #elif mods=='ocean_types':
    #    suffixes=['SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt','SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_oligotrophic.txt', 'SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_upwelling.txt']
    #    titles=['SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic','SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_oligotrophic', 'SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_upwelling']
    
    files_algae=['model_data/Tot_algae_'+suff for suff in suffixes]
    files_virus=['model_data/Virus_'+suff for suff in suffixes]
    files_zoop=['model_data/Zoop_'+suff for suff in suffixes]
    files_ZV=['model_data/ZV_ratio_'+suff for suff in suffixes]
    files_S=['model_data/Susceptible_'+suff for suff in suffixes]
    files_I=['model_data/Infected_'+suff for suff in suffixes]
    files_NPP=['model_data/NPP_'+suff for suff in suffixes]
    if mods=='main_models':
        files_R=['model_data/Resistant_'+suffixes[i] for i in range(1,3)]

    def read_data(files, mi_h):
        datas=[]
        mxs=[]
        mis=[]
        for i in range(len(files)):
            f=files[i]
            u=load_matrix(f, ' ')
            u=np.array(u)
            mx=np.nanmax(u)
            mi=np.nanmin(u)
            datas.append(u)
            mxs.append(mx)
            mis.append(mi)
        mx=max(mxs)
        mi=min(mis)

        mi_f=max([mi_h, mi])

        return datas, mx, mi_f

    # load data
    datas_a, mx_a, mi_a=read_data(files_algae, 5)
    datas_v, mx_v, mi_v=read_data(files_virus, 5)
    datas_z, mx_z, mi_z=read_data(files_zoop, 3)
    datas_zv, mx_zv, mi_zv=read_data(files_ZV, -4)
    datas_s, mx_s, mi_s=read_data(files_S, -1000)
    datas_i, mx_i, mi_i=read_data(files_I, 0)
    datas_npp, mx_npp, mi_npp=read_data(files_NPP, -1.5)
    if mods=='main_models':
        datas_r, mx_r, mi_r=read_data(files_R, -1000)

    phis=list(np.logspace(-12, -8, (12-8)*9+1))
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]

    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]

    ylb='Latent period'
    # pdf for tracers concentrations, ZV ratio and NPP
    pp=PdfPages('Scaled_PVZ_'+mods+'.pdf')
    for j,d in enumerate(datas_a):
        mi=mi_a
        mx=mx_a
        plot_with_scale_bis(d, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S+I (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_v):
        mi=mi_v
        mx=mx_v
        plot_with_scale_bis(d, 'Blues', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'V (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_z):
        mi=mi_z
        mx=mx_z
        plot_with_scale_bis(d, 'Reds', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Z (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_zv):
        mi=mi_zv
        mx=mx_zv
        amx=max([abs(mi), abs(mx)])
        cmap=matplotlib.colormaps['bwr']
        rcmap = cmap.reversed()
        plot_with_scale_bis(d, rcmap, -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky , 'ZV ratio (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_s):
        mi=mi_s
        mx=mx_s
        plot_with_scale_bis(d, 'Greens', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'S (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_i):
        mi=mi_i
        mx=mx_i
        plot_with_scale_bis(d, 'Oranges', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'I (log10(ind.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_npp):
        mi=mi_npp
        mx=mx_npp
        cmap = plt.cm.RdYlGn
        plot_with_scale_bis(d, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'NPP (log10(umolC.d-1.L-1)) '+titles[j], yl=ylb)
        pp.savefig()
    if mods=='main_models':
        for j in range(1,3):
            d=datas_r[j-1]
            mi=mi_r
            mx=mx_r
            plot_with_scale_bis(d, 'Purples', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'R (log10(ind.L-1)) '+titles[j], yl=ylb)
            pp.savefig()

    pp.close()


    
    # Fourier analysis
    suffixes=['SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.txt','SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt','SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']
    titles=['SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic','SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic','SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic', 'SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic']

    files_FFT=['model_data/FFT_modulo_'+suff for suff in suffixes]
    files_period=['model_data/Period_'+suff for suff in suffixes]

    datas_fft, mx_fft, mi_fft=read_data(files_FFT, -1000)
    datas_per, mx_per, mi_per=read_data(files_period, -1000)

    mi_pers=[]
    mx_pers=[]
    for i, dat in enumerate(datas_per):
        ff=datas_fft[i]
        dat[ff<-2]=np.nan
        mi_pers.append(np.nanmin(dat))
        mx_pers.append(np.nanmax(dat))
        datas_per[i]=dat

    mx_per, mi_per=np.nanmax(np.array(mx_pers)), np.nanmin(np.array(mi_pers))

    # Fourier pdf
    pp=PdfPages('Scaled_PVZ_FFT_and_Period.pdf')
    for j,d in enumerate(datas_fft):
        mi=mi_fft
        mx=mx_fft
        plot_with_scale_bis(d, 'inferno', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Modulo FFT '+titles[j], yl=ylb)
        pp.savefig()
    for j,d in enumerate(datas_per):
        mi=mi_per
        mx=mx_per
        plot_with_scale_bis(d, 'plasma', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , 'Period '+titles[j], yl=ylb)
        pp.savefig()
    pp.close()
        
    
    files_coexistence_SIVZ= ['coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Fluctuation_free_growth_rate.txt', 'coexistence_analysis_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic_Relative_non_linearity_I.txt']
    data_coex, mx_coex, mi_coex=read_data(files_coexistence_SIVZ, -100000000000000000000000000)

    state_reached=load_matrix('model_data/state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.txt', ' ')
    state_reached=np.array(state_reached)

    coex_mat=copy.deepcopy(state_reached)
    coex_mat[state_reached!=4]=0
    coex_mat[state_reached==4]=1
    coex_mat=np.transpose(coex_mat)
    coex_mat=np.flipud(coex_mat)

    # Coexistence analysis
    pp=PdfPages('Scaled_SIVZ_coexistence.pdf')
    for j,d in enumerate(data_coex):
        dat=abs(d)
        dat[dat<1e-1]=1e-1
        cmap = matplotlib.colormaps['PuOr_r']
        colors0 = cmap(np.linspace(0, 1, 200))
        if j==0:
            colors0=colors0[100:]
        else:
            colors0=colors0[100:0:-1]
        cmap0=matplotlib.colors.ListedColormap(colors0)
        mi=mt.log10(np.nanmin(dat))
        mx=mt.log10(np.nanmax(dat))
        bounds=np.linspace(mi, mx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        axi=plot_with_scale_bis(np.log10(dat), cmap0, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , files_coexistence_SIVZ[j], yl=ylb)
        draw_borders(coex_mat, 1.5, 0.02, 'black', axi)
        pp.savefig()
    for j,d in enumerate(data_coex):
        dat=abs(d)
        dat[dat>100]=100#
        cmap = matplotlib.colormaps['PuOr_r']
        colors0 = cmap(np.linspace(0, 1, 200))
        if j==0:
            colors0=colors0[100:]
        else:
            colors0=colors0[100:0:-1]
        cmap0=matplotlib.colors.ListedColormap(colors0)
        mi=mt.log10(np.nanmin(dat))
        mx=mt.log10(np.nanmax(dat))
        bounds=np.linspace(mi, mx, 100)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
        axi=plot_with_scale_bis(np.log10(dat), cmap0, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky , files_coexistence_SIVZ[j], yl=ylb)
        draw_borders(coex_mat, 1.5, 0.02, 'black', axi)
        pp.savefig()
    pp.close()


    files_distances=['model_data/distance_A_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt','model_data/distance_A_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/distance_A_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt', 'model_data/distance_V_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt','model_data/distance_V_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/distance_V_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt', 'model_data/distance_Z_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt','model_data/distance_Z_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/distance_Z_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt']
    A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc=concentration_ranges(3, 'mesotrophic')

    data_dist, mx_dist, mi_dist=read_data(files_distances, -1000000)

    data_distances=[]
    mx_distances=[]
    mi_distances=[]
    for j,d in enumerate(data_dist):
        if j<3:
            scaled_d=d*100/target_conc[0]
        elif j>=3 and j<6:
            scaled_d=d*100/target_conc[1]
        elif j>=6 and j<9:
            scaled_d=d*100/target_conc[2]
        data_distances.append(scaled_d)
        mx_distances.append(np.nanmax(scaled_d))
        mi_distances.append(np.nanmin(scaled_d))

    mx_distance=max(mx_distances)
    mi_distance=min(mi_distances)
    amx=max([abs(mi_distance),abs(mx_distance)])
    amx=mt.ceil(amx/ 10) * 10
    # distances to specififc targets (A, Z and V) (A stands for Algae/Phytoplankton)
    pp=PdfPages('Scaled_SIVZ_distances.pdf')
    for j,d in enumerate(data_distances):
        plot_with_scale_bis(d, 'coolwarm', -amx, amx, atickx, aticky, alabel_tickx, alabel_ticky, files_distances[j], yl=ylb)
        pp.savefig()
    pp.close()


    files_distances=['model_data/distances_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt','model_data/distances_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/distances_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt' ]
    data_dist, mx_dist, mi_dist=read_data(files_distances, -1000000)

    mx_distance=mt.ceil(mx_dist/ 10) * 10
    mi_distance=0
    amx=mx_distance
    pp=PdfPages('Scaled_SIVZ_distances_tot.pdf')
    for j,d in enumerate(data_dist):
        plot_with_scale_bis(d, 'inferno', 0, amx, atickx, aticky, alabel_tickx, alabel_ticky, files_distances[j], yl=ylb)
        pp.savefig()
    pp.close()

