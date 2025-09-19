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
def make_state_diagram_simulated(final_Surv, f):
    if 'SVZ' in f or 'SVRZ' in f and 'mesotrophic_phir.txt' not in f and 'mesotrophic_epsr.txt' not in f:
        final_Surv = final_Surv[np.newaxis, :]
        final_Surv=np.transpose(final_Surv)
    if 'SIVZ' in f or 'SVZ' in f:
        bounds=[0,1,2,3,4]
        norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)
        mi=0
        mx=4
        idis=[0,2,6,4,9]
    else:
        mi=0
        mx=10
        bounds=[0,1,2,3,4,5,6,7,8,9,10]
        norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)

    if 'SVZ' in f or 'SVRZ' in f:
        if 'mesotrophic_phir.txt' in f or 'mesotrophic_epsr.txt' in f:
            atickx=[i*9 if i>0 else 0 for i in range(5) ]
            aticky=[i*9+9 if i!=0 else i*9+9 for i in range(3,-1,-1)]
            alabel_tickx= [phis[i] for i in atickx]
            atickyr=[i*9+8 if i!=0 else i*9+9 for i in range(0,4)]
            alabel_ticky= [phi_over_phir[i-9] for i in atickyr]
            aticky[0]=aticky[0]-1
            if 'mesotrophic_phir.txt' in f:
                ylb='phi/phir'
            elif 'mesotrophic_epsr.txt' in f:
                ylb='eps/epsr'
        else:
            atickx=[i*9 if i>0 else 0 for i in range(5)]
            aticky=[]
            alabel_tickx= [phis[i] for i in atickx]
            atickyr=Reverse(aticky)
            alabel_ticky= []
            ylb=''
    else:
        if 'epsr_mesotrophic.txt' in f or 'phir_mesotrophic.txt' in f:
            atickx=[i*9 if i>0 else 0 for i in range(5) ]
            aticky=[i*9+9 if i!=0 else i*9+9 for i in range(3,-1,-1)]
            alabel_tickx= [phis[i] for i in atickx]
            atickyr=[i*9+8 if i!=0 else i*9+9 for i in range(0,4)]
            alabel_ticky= [phi_over_phir[i-9] for i in atickyr]
            aticky[0]=aticky[0]-1
            if 'epsr_mesotrophic.txt' in f:
                ylb='eps/epsr'
            elif 'phir_mesotrophic.txt' in f:
                ylb='phi/phir'
        else:
            atickx=[i*9 if i>0 else 0 for i in range(5) ]
            aticky=[i*9+8 for i in range(2,-1,-1) ]
            alabel_tickx= [phis[i] for i in atickx]
            atickyr=Reverse(aticky)
            alabel_ticky= [lps[i-8] for i in atickyr]
            ylb='Latent period'
    
    n_state_R=11
    mx_R=n_state_R-1
    colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx_R+2)])
    cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
    colors.append(cn)
    colors[0]=(0,0,0)
    colors[10]=(0.8039, 0.5216, 0.2471)
    if 'SIVZ' in f or 'SVZ' in f:
        colors=[colors[idi] for idi in idis]
    colors=tuple(colors)
    cmap_R= matplotlib.colors.ListedColormap(colors)
    axi=plot_with_scale(final_Surv,cmap_R,mi,mx, atickx, aticky, alabel_tickx, alabel_ticky, 'State reached', norm=norm, yl=ylb)
    #print(final_Surv.shape)
    coex_mat=np.zeros( final_Surv.shape )
    #print(final_Surv)
    if 'SIVZ' in f or 'SVZ' in f:
        coex_mat[final_Surv==4]=1
    else:
        coex_mat[final_Surv>7]=1
    coex_mat=np.transpose(coex_mat)
    coex_mat=np.flipud(coex_mat)
    #print(coex_mat)
    if 'SVZ' in f or 'SVRZ' in f and 'mesotrophic_phir.txt' not in f and 'mesotrophic_epsr.txt' not in f:
        draw_borders(coex_mat, 1.5, 0.02, 'white', axi)
    else:
        draw_borders(coex_mat, 4.5, 0.3, 'white', axi)
    pp.savefig()

def make_state_diagram_th(predicted_state_bis,alter_state, u_state,limit_val_mat, f):
    atickx=[i*9 if i>0 else 0 for i in range(5) ]
    aticky=[i*9+8 for i in range(2,-1,-1) ]
    alabel_tickx= [phis[i] for i in atickx]
    atickyr=Reverse(aticky)
    alabel_ticky= [lps[i-8] for i in atickyr]
    ylb='Latent period'

    if 'SIVZ' in f or 'SVZ' in f:
        bounds=[0,1,2,3,4]
        norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)
        mi=0
        mx=4
        idis=[0,2,6,4,9]
    else:
        mi=0
        mx=10
        bounds=[0,1,2,3,4,5,6,7,8,9,10]
        norm = matplotlib.colors.BoundaryNorm(bounds, len(bounds)-1)

    n_state_R=11
    mx_R=n_state_R-1
    colors = list(matplotlib.colormaps.get_cmap('tab20').colors[2:(mx_R+2)])
    cn= matplotlib.colormaps.get_cmap('tab20').colors[0]
    colors.append(cn)
    colors[0]=(0,0,0)
    colors[10]=(0.8039, 0.5216, 0.2471)
    if 'SIVZ' in f or 'SVZ' in f:
        colors=[colors[idi] for idi in idis]
    colors=tuple(colors)
    cmap_R= matplotlib.colors.ListedColormap(colors)
        
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
    if 'SIVZ' in f or 'SVZ' in f:
        coex_mat[predicted_state_bis==4]=1
    else:
        coex_mat[predicted_state_bis>7]=1
    coex_mat=np.transpose(coex_mat)
    coex_mat=np.flipud(coex_mat)
    draw_borders(coex_mat, 4.5, 0.3, 'white', axi)
    pp.savefig()



if __name__ == '__main__':
    files_suffixes= ['state_reached_SVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.txt', 'state_reached_SVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'state_reached_SVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_m2-1800.0_mesotrophic.txt', 'state_reached_SVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_no-dz2_mesotrophic.txt', 'state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_no-dz2_mesotrophic.txt', 'state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-5.0.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-5.0.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-5.0.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800.0_mesotrophic_phir-5.0.txt' , 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-5.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt','state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir-fullres.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir-fullres.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir-fullres.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800.0_mesotrophic_phir-fullres.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_lp_phi-ratio-0.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-0.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-0.0_mesotrophic.txt']
    files=['model_data/'+f for f in files_suffixes]

    #ylb='Latent period'
    phis=list(np.logspace(-12, -8, (12-8)*9+1))
    lps=list(np.logspace(-1, 2, (2+1)*9+1))
    lps=lps[0:len(lps)-1]

    phi_over_phir=list(np.logspace(0, 4, 4*9+1))
    phi_over_phir=phi_over_phir[0:len(phi_over_phir)-1]
    phi_over_phir.append(0)
    del phi_over_phir[0]

    pp = PdfPages('State_diagrams_main_models.pdf')
    for f in files:
        print(f)
        final_Surv, mxx, mii=read_data([f], -10)
        #print(final_Surv)
        
        final_Surv=np.array(final_Surv)
        #if 'SVZ' in f:
        #final_Surv = final_Surv.squeeze(-1)  
        final_Surv = np.squeeze(final_Surv)
        make_state_diagram_simulated(final_Surv, f) 
    pp.close()

    files_suffixes_res=['state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_phir.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_phir.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_phir.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800.0_mesotrophic_phir.txt','state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_phir_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_phir_mesotrophic.txt','state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_phir_mesotrophic.txt','state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_phir_mesotrophic.txt' ,'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_mesotrophic_epsr.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_mesotrophic_epsr.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_mesotrophic_epsr.txt', 'state_reached_SVRZ_Prochlorochoccus_BS15.0_LOI0_GR-R0.8_m2-1800.0_mesotrophic_epsr.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_epsr_mesotrophic.txt', 'state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_epsr_mesotrophic.txt','state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_no-dz2_m2-1800.0_epsr_mesotrophic.txt','state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_epsr_mesotrophic.txt']
    files_res=['model_data/'+f for f in files_suffixes_res]

    pp = PdfPages('State_diagrams_resistance_models.pdf')
    for f in files_res:
        print(f)
        final_Surv, mxx, mii=read_data([f], -10)
        #print(final_Surv)

        final_Surv=np.array(final_Surv)
        #if 'SVZ' in f:
        #final_Surv = final_Surv.squeeze(-1)
        
        final_Surv = np.squeeze(final_Surv)
        #print(final_Surv.shape)
        make_state_diagram_simulated(final_Surv, f)
    pp.close()

    files_PT=['model_data/state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Synechococcus_BS30.0_LOI0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Synechococcus_BS30.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Eukaryote_BS180.0_LOI0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Eukaryote_BS180.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Diatom_BS270.0_LOI0_mesotrophic.txt', 'model_data/state_reached_SIVZ_Diatom_BS270.0_LOI0_m2-1800.0_mesotrophic.txt']

    pp = PdfPages('State_diagrams_phytoplankton_types_models.pdf')
    for f in files_PT:
        print(f)
        final_Surv, mxx, mii=read_data([f], -10)
        #print(final_Surv)

        final_Surv=np.array(final_Surv)
        #if 'SVZ' in f:
        #final_Surv = final_Surv.squeeze(-1)
        if 'Diatom' in f:
            phis=list(np.logspace(-11, -7, (12-8)*9+1))
        final_Surv = np.squeeze(final_Surv)
        make_state_diagram_simulated(final_Surv, f)
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    pp.close()

    files_PT_res=['model_data/state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt','model_data/state_reached_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/state_reached_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']

    pp = PdfPages('State_diagrams_phytoplankton_types_resistance_models.pdf')
    for f in files_PT_res:
        print(f)
        final_Surv, mxx, mii=read_data([f], -10)

        final_Surv=np.array(final_Surv)
        if 'Diatom' in f:
            phis=list(np.logspace(-11, -7, (12-8)*9+1))
        final_Surv = np.squeeze(final_Surv)
        make_state_diagram_simulated(final_Surv, f)
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    pp.close()

    files_PT_th=['model_data/predicted_state_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Synechococcus_BS30.0_LOI0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Synechococcus_BS30.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Eukaryote_BS180.0_LOI0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Eukaryote_BS180.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Diatom_BS270.0_LOI0_mesotrophic.txt', 'model_data/predicted_state_th_SIVZ_Diatom_BS270.0_LOI0_m2-1800.0_mesotrophic.txt']
    files_PT_th_alt=['model_data/alternate_stable_state_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Synechococcus_BS30.0_LOI0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Synechococcus_BS30.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Eukaryote_BS180.0_LOI0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Eukaryote_BS180.0_LOI0_m2-1800.0_mesotrophic.txt','model_data/alternate_stable_state_th_SIVZ_Diatom_BS270.0_LOI0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVZ_Diatom_BS270.0_LOI0_m2-1800.0_mesotrophic.txt']
    files_PT_th_unst=['model_data/unstable_equilibria_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Synechococcus_BS30.0_LOI0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Synechococcus_BS30.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Eukaryote_BS180.0_LOI0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Eukaryote_BS180.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Diatom_BS270.0_LOI0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVZ_Diatom_BS270.0_LOI0_m2-1800.0_mesotrophic.txt']
    files_PT_th_lim=['model_data/limit_val_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Prochlorochoccus_BS15.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Synechococcus_BS30.0_LOI0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Synechococcus_BS30.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Eukaryote_BS180.0_LOI0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Eukaryote_BS180.0_LOI0_m2-1800.0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Diatom_BS270.0_LOI0_mesotrophic.txt', 'model_data/limit_val_th_SIVZ_Diatom_BS270.0_LOI0_m2-1800.0_mesotrophic.txt']
    
    pp = PdfPages('State_diagrams_phytoplankton_types_th.pdf')
    for i, f in enumerate(files_PT_th):
        pred, mxx, mii=read_data([f], -10)
        alt, mxx, mii=read_data([files_PT_th_alt[i]], -10)
        unst, mxx, mii=read_data([files_PT_th_unst[i]], -10)
        lim, mxx, mii=read_data([files_PT_th_lim[i]], -10)

        pred=np.array(pred)
        alt=np.array(alt)
        unst=np.array(unst)
        lim=np.array(lim)
        if 'Diatom' in f:
            phis=list(np.logspace(-11, -7, (12-8)*9+1))
        pred=np.squeeze(pred)
        alt=np.squeeze(alt)
        unst=np.squeeze(unst)
        lim=np.squeeze(lim)
        make_state_diagram_th(pred,alt, unst, lim, f)
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    pp.close()


    files_PT_res_th=['model_data/predicted_state_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt',  'model_data/predicted_state_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/predicted_state_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']
    files_PT_res_th_alt=['model_data/alternate_stable_state_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/alternate_stable_state_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']
    files_PT_res_th_unst=['model_data/unstable_equilibria_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt','model_data/unstable_equilibria_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/unstable_equilibria_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']
    files_PT_res_th_lim=['model_data/limit_val_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Prochlorochoccus_LP0.37_BS15.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Synechococcus_LP0.37_BS30.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Eukaryote_LP0.37_BS180.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_lp_phi-ratio-5.0_mesotrophic.txt', 'model_data/limit_val_th_SIVRZ_Diatom_LP0.95_BS270.0_LOI0_GR-R0.8_m2-1800.0_lp_phi-ratio-5.0_mesotrophic.txt']
    pp = PdfPages('State_diagrams_phytoplankton_types_resistance_th.pdf')
    for i, f in enumerate(files_PT_res_th):
        pred, mxx, mii=read_data([f], -10)
        alt, mxx, mii=read_data([files_PT_res_th_alt[i]], -10)
        unst, mxx, mii=read_data([files_PT_res_th_unst[i]], -10)
        lim, mxx, mii=read_data([files_PT_res_th_lim[i]], -10)

        pred=np.array(pred)
        alt=np.array(alt)
        unst=np.array(unst)
        lim=np.array(lim)
        if 'Diatom' in f:
            phis=list(np.logspace(-11, -7, (12-8)*9+1))
        pred=np.squeeze(pred)
        alt=np.squeeze(alt)
        unst=np.squeeze(unst)
        lim=np.squeeze(lim)
        make_state_diagram_th(pred,alt, unst, lim, f)
        phis=list(np.logspace(-12, -8, (12-8)*9+1))
    pp.close()

