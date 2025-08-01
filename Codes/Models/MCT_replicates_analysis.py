import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import sys
import copy
from numpy import inf
from scipy.signal import find_peaks
from generic_functions import *


if __name__ == '__main__':
    data_replicates_chesson=[]
    data_replicates_invasion_growth=[]
    data_replicates_invasion_V=[]
    data_replicates_invasion_Z=[]
    data_replicates_invasion_Z_ffree=[]
    data_replicates_invasion_Z_S=[]
    data_replicates_invasion_V_ffree=[]
    for i in range(10):
        dat_chesson=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_chesson_binary.txt', ' ')
        dat_invasion_V=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_invasion_growth_rate_V_inv.txt', ' ')
        dat_invasion_Z=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_invasion_growth_rate_Z_inv.txt', ' ' )
        dat_invasion_growth=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_chesson_growth.txt', ' ')
        dat_invasion_Z_ffree=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_fluctuation_free_growth_rate_Z_inv.txt', ' ')
        dat_invasion_Z_S=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_relative_non_linearity_in_S_Z_inv.txt', ' ')
        dat_invasion_V_ffree=load_matrix('model_data/SIVZ_MCT_Prochlorochoccus_BS15.0_LOI0_replicate_'+str(i)+'_no-dz2_no-dv2_mesotrophic_fluctuation_free_growth_rate_V_inv.txt', ' ')
    
        data_replicates_chesson.append(dat_chesson)
        data_replicates_invasion_growth.append(dat_invasion_growth)
        data_replicates_invasion_V.append(dat_invasion_V)
        data_replicates_invasion_Z.append(dat_invasion_Z)
        data_replicates_invasion_Z_ffree.append(dat_invasion_Z_ffree)
        data_replicates_invasion_Z_S.append(dat_invasion_Z_S)
        data_replicates_invasion_V_ffree.append(dat_invasion_V_ffree)

    state_reached=load_matrix('model_data/state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.txt', ' ')
    state_reached=np.array(state_reached)

    def percentage_of_ones(array_list):
        stacked = np.stack(array_list, axis=0)
        percent_ones = 100 * stacked.mean(axis=0)

        return percent_ones

    def mean_list_arrays(array_list):
        mean_array = np.mean(array_list, axis=0)
        return mean_array
    
    def percentage_positive(array_list):
        stacked = np.stack(array_list, axis=0)  # shape: (n_arrays, ...)
        positive_mask = stacked > 0
        percent_positive = 100 * positive_mask.mean(axis=0)
        return percent_positive

    percent_ones=percentage_of_ones(data_replicates_chesson)
    all_invasion_positive = (percent_ones==100).astype(float)

    preds_chesson1=all_invasion_positive[state_reached==4.0]
    true_values1=[1 for i in range(len(preds_chesson1))]
    preds_chesson1=np.array(preds_chesson1)
    true_values1=np.array(true_values1)
    true_positives1 = np.sum((true_values1 == 1) & (preds_chesson1 == 1))
    actual_positives1 = np.sum(true_values1 == 1)
    sensitivity1=true_positives1 / actual_positives1
    print('Sensitivity of Chesson criterion with SIV collapse region')
    print(sensitivity1)

    all_invasion_positive[np.isnan(percent_ones)]=np.nan
    mean_invasion_growth=mean_list_arrays(data_replicates_invasion_growth)
    mask = (np.isclose(state_reached, 4.0)) & (~np.isnan(all_invasion_positive))
    preds_chesson2=all_invasion_positive[mask]
    preds_chesson2=np.array(preds_chesson2)
    true_values2=[1 for i in range(len(preds_chesson2))]
    true_values2=np.array(true_values2)
    true_positives2 = np.sum((true_values2 == 1) & (preds_chesson2 == 1))
    actual_positives2 = np.sum(true_values2 == 1)
    sensitivity2=true_positives2 / actual_positives2
    print('Sensitivity of Chesson criterion without SIV collapse region')
    print(sensitivity2)

    mean_invasion_V=mean_list_arrays(data_replicates_invasion_V)
    mean_invasion_Z=mean_list_arrays(data_replicates_invasion_Z)
    percent_supp0_V=percentage_positive(data_replicates_invasion_V)
    percent_supp0_Z=percentage_positive(data_replicates_invasion_Z)
    all_supp0_V=(percent_supp0_V==100).astype(int)
    all_supp0_Z=(percent_supp0_Z==100).astype(int)
    one_supp0_V=(percent_supp0_V>0).astype(int)
    one_supp0_Z=(percent_supp0_Z>0).astype(int)

    mean_invasion_Z_ffree=mean_list_arrays(data_replicates_invasion_Z_ffree)
    mean_invasion_Z_S=mean_list_arrays(data_replicates_invasion_Z_S)
    mean_invasion_V_ffree=mean_list_arrays(data_replicates_invasion_V_ffree)
    
    percent_supp0_Z_ffree=percentage_positive(data_replicates_invasion_Z_ffree)
    percent_supp0_Z_S=percentage_positive(data_replicates_invasion_Z_S)
    percent_supp0_V_ffree=percentage_positive(data_replicates_invasion_V_ffree)

    all_supp0_Z_ffree=(percent_supp0_Z_ffree==100).astype(int)
    all_supp0_Z_S=(percent_supp0_Z_S==100).astype(int)
    all_supp0_V_ffree=(percent_supp0_V_ffree==100).astype(int)
    one_supp0_Z_ffree=(percent_supp0_Z_ffree>0).astype(int)
    one_supp0_Z_S=(percent_supp0_Z_S>0).astype(int)
    one_supp0_V_ffree=(percent_supp0_V_ffree>0).astype(int)
    
    # adsorption rate parameter space
    phis=list(np.logspace(-12, -8, (12-8)*9+1))

    #latent period parameter space
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]
    
    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]
    ylb='Latent period'

    pp = PdfPages('SIVZ_model_phi_latent_period_replicates_consensus_MCT.pdf')

    mi=0
    mx=100
    plot_with_scale_bis(percent_ones, 'nipy_spectral', mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, '% of valid chesson accross replicates', yl=ylb)
    pp.savefig()


    mi=0
    mx=1
    cmap = matplotlib.colors.ListedColormap(['white', 'green'])
    plot_with_scale(all_invasion_positive, cmap, mi, mx, atickx, aticky, alabel_tickx, alabel_ticky, 'Invasibility (valid Chesson across all replicates)')
    pp.savefig()

    cmap = matplotlib.colormaps['PuOr_r']
    colors0 = cmap(np.linspace(0, 1, 100))
    cmap0=matplotlib.colors.ListedColormap(colors0)

    mx=np.nanmax(np.array(mean_invasion_V))
    mi=np.nanmin(np.array(mean_invasion_V))
    amx=np.nanmax(np.array([abs(mi), abs(mx)]))
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=-amx
    ma=amx
    axi=plot_with_scale(mean_invasion_V, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion growth rate (V) + all supp0', norm=norm)
    all_supp0_V=np.transpose(all_supp0_V)
    all_supp0_V=np.flipud(all_supp0_V)
    draw_borders(all_supp0_V, 4.5, 0.2, 'green', axi)
    pp.savefig()

    axi=plot_with_scale(mean_invasion_V, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion growth rate (V) + all supp0 + at least 1 supp0', norm=norm)
    draw_borders(all_supp0_V, 4.5, 0.2, 'green', axi)
    one_supp0_V=np.transpose(one_supp0_V)
    one_supp0_V=np.flipud(one_supp0_V)
    draw_borders(one_supp0_V, 4.5, 0.2, 'blue', axi)
    pp.savefig()

    mx=np.nanmax(np.array(mean_invasion_Z))
    mi=np.nanmin(np.array(mean_invasion_Z))
    amx=np.nanmax(np.array([abs(mi), abs(mx)]))
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=-amx
    ma=amx
    axi=plot_with_scale(mean_invasion_Z, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion growth rate (Z) + all supp0', norm=norm)
    all_supp0_Z=np.transpose(all_supp0_Z)
    all_supp0_Z=np.flipud(all_supp0_Z)
    draw_borders(all_supp0_Z, 4.5, 0.2, 'green', axi)
    pp.savefig()

    axi=plot_with_scale(mean_invasion_Z, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion growth rate (Z) + all supp0 + at least 1 supp0', norm=norm)
    draw_borders(all_supp0_Z, 4.5, 0.2, 'green', axi)
    one_supp0_Z=np.transpose(one_supp0_Z)
    one_supp0_Z=np.flipud(one_supp0_Z)
    draw_borders(one_supp0_Z, 4.5, 0.2, 'blue', axi)
    pp.savefig()

    mx=np.nanmax(np.array(mean_invasion_Z_ffree))
    mi=np.nanmin(np.array(mean_invasion_Z_ffree))
    amx=np.nanmax(np.array([abs(mi), abs(mx)]))
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=-amx
    ma=amx
    axi=plot_with_scale(mean_invasion_Z_ffree, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion fluctuation free growth rate (Z inv) + all supp0', norm=norm)
    all_supp0_Z_ffree=np.transpose(all_supp0_Z_ffree)
    all_supp0_Z_ffree=np.flipud(all_supp0_Z_ffree)
    draw_borders(all_supp0_Z_ffree, 4.5, 0.2, 'green', axi)
    pp.savefig()

    axi=plot_with_scale(mean_invasion_Z_ffree, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion fluctuation free growth growth rate (Z inv) + all supp0 + at least 1 supp0', norm=norm)
    draw_borders(all_supp0_Z_ffree, 4.5, 0.2, 'green', axi)
    one_supp0_Z_ffree=np.transpose(one_supp0_Z_ffree)
    one_supp0_Z_ffree=np.flipud(one_supp0_Z_ffree)
    draw_borders(one_supp0_Z_ffree, 4.5, 0.2, 'blue', axi)
    pp.savefig()

    mx=np.nanmax(np.array(mean_invasion_Z_S))
    mi=np.nanmin(np.array(mean_invasion_Z_S))
    amx=np.nanmax(np.array([abs(mi), abs(mx)]))
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=-amx
    ma=amx
    axi=plot_with_scale(mean_invasion_Z_S, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average relative non linearity in S (Z inv) + all supp0', norm=norm)
    all_supp0_Z_S=np.transpose(all_supp0_Z_S)
    all_supp0_Z_S=np.flipud(all_supp0_Z_S)
    draw_borders(all_supp0_Z_S, 4.5, 0.2, 'green', axi)
    pp.savefig()

    axi=plot_with_scale(mean_invasion_Z_S, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average relative non linearity in S (Z inv) + all supp0 + at least 1 supp0', norm=norm)
    draw_borders(all_supp0_Z_S, 4.5, 0.2, 'green', axi)
    one_supp0_Z_S=np.transpose(one_supp0_Z_S)
    one_supp0_Z_S=np.flipud(one_supp0_Z_S)
    draw_borders(one_supp0_Z_S, 4.5, 0.2, 'blue', axi)
    pp.savefig()

    mx=np.nanmax(np.array(mean_invasion_V_ffree))
    mi=np.nanmin(np.array(mean_invasion_V_ffree))
    amx=np.nanmax(np.array([abs(mi), abs(mx)]))
    bounds=np.linspace(-amx, amx, 100)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap0.N-1)
    mi=-amx
    ma=amx
    axi=plot_with_scale(mean_invasion_V_ffree, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion fluctuation free growth rate (V inv) + all supp0', norm=norm)
    all_supp0_V_ffree=np.transpose(all_supp0_V_ffree)
    all_supp0_V_ffree=np.flipud(all_supp0_V_ffree)
    draw_borders(all_supp0_V_ffree, 4.5, 0.2, 'green', axi)
    pp.savefig()

    axi=plot_with_scale(mean_invasion_V_ffree, cmap0,mi,ma, atickx, aticky, alabel_tickx, alabel_ticky, 'Average invasion fluctuation free growth growth rate (Z inv) + all supp0 + at least 1 supp0', norm=norm)
    draw_borders(all_supp0_V_ffree, 4.5, 0.2, 'green', axi)
    one_supp0_V_ffree=np.transpose(one_supp0_V_ffree)
    one_supp0_V_ffree=np.flipud(one_supp0_V_ffree)
    draw_borders(one_supp0_V_ffree, 4.5, 0.2, 'blue', axi)
    pp.savefig()
 
    
    pp.close()


