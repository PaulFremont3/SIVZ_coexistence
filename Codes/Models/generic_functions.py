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
from scipy.signal import find_peaks

# function defining realistic concentrations
def concentration_ranges(indice, otype):
    if indice in [2,3] and otype in ['upwelling', '0']:
        A_cond_low=1e8
        A_cond_high=1e10
        V_cond_low=1e8
        V_cond_high=1e11
        Z_cond_low=1e3
        Z_cond_high=1e6
        I_cond_high=100
        I_cond_low=5
        perc_cond_high=100
        perc_cond_low=0
    if indice in [0,1] and otype in ['upwelling', '0']:
        A_cond_low=1e6
        A_cond_high=1e8
        V_cond_low=1e8
        V_cond_high=1e10
        Z_cond_low=1e3
        Z_cond_high=1e6
        I_cond_high=100
        I_cond_low=5
        perc_cond_high=100
        perc_cond_low=0
    if indice in [2,3] and otype in ['oligotrophic', 'mesotrophic']:
        A_cond_low=1e7
        A_cond_high=1e9
        V_cond_low=1e7
        V_cond_high=1e10
        Z_cond_low=1e3
        Z_cond_high=1e6
        I_cond_high=10
        if indice==0:
            I_cond_high=50
        I_cond_low=0.5
        perc_cond_high=50
        perc_cond_low=0
    if indice in [0,1] and otype in ['oligotrophic', 'mesotrophic']:
        if indice==0:
            A_cond_low=1e5
            A_cond_high=5e7
        elif indice==1:
            A_cond_low=1e6
            A_cond_high=1e8
        if indice==1:
            V_cond_low=1e7
            V_cond_high=1e9
        elif indice==0:
            V_cond_low=1e8
            V_cond_high=1e10
        Z_cond_low=1e3
        Z_cond_high=1e6
        I_cond_high=10
        if indice==0:
            I_cond_high=50
        I_cond_low=0.5
        perc_cond_high=50
        perc_cond_low=0

    target_conc_u, target_conc_m, target_conc_o=target_concentrations(indice)
    if otype=='oligotrophic':
        target_conc=target_conc_o
    elif otype=='mesotrophic':
        target_conc=target_conc_m
    elif otype=='upwelling':
        target_conc=target_conc_u
    return A_cond_low, A_cond_high, V_cond_low, V_cond_high, Z_cond_low, Z_cond_high, I_cond_high, I_cond_low, perc_cond_high, perc_cond_low, target_conc

# target concentrations
def target_concentrations(indice):
    if indice ==3:
        target_conc_u=[4e8, 2e9, 2e5]
        target_conc_m=[2e8, 1e9, 1e5]
        target_conc_o=[1e8, 6e8, 1e4]
    if indice==2:
        target_conc_u=[4e7, 2e9, 2e5]
        target_conc_m=[5e7, 1e9, 1e5]
        target_conc_o=[2e7, 6e8, 1e4]
    if indice==0:
        target_conc_u=[4e6, 2e9, 1e5]
        target_conc_m=[2e6, 1e9, 5e4]
        target_conc_o=[5e5, 2e8, 5e3]
    if indice==1:
        target_conc_u=[4e7, 4e8, 1e5]
        target_conc_m=[2e7, 2e8, 5e4]
        target_conc_o=[1e7, 1e8, 5e3]
    return target_conc_u, target_conc_m, target_conc_o

# to define stability based on jacobian eigenvalues
def check_eigs(eigs, eig_vecs, ast, ost):
    ast=1
    ost=1
    for ei in eigs:
        ei_r=ei.real
        ei_i=ei.imag
        if ei_r>0:
            ast=0
        if ei_i!=0:
            ost=0
    if ast==0:
        fst=0
    elif ast==1 and ost==0:
        fst=1
    elif ast==1 and ost==1:
        fst=2
    return ast, ost, fst

# Reversing a list
def Reverse(lst):
    new_lst = lst[::-1]
    return new_lst

# plot functions of matrices
def plot_without_scale(matrix, tickx, ticky, label_tickx, label_ticky, title, yl):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.imshow(np.flipud(np.transpose(matrix)))
    ax.set_xticks(tickx)
    ax.set_xticklabels(label_tickx, rotation='vertical')
    ax.set_yticks(ticky)
    ax.set_yticklabels(label_ticky)
    plt.xlabel('Phi')
    plt.ylabel(yl)
    plt.title(title)

def plot_with_scale(matrix, colorcode, mi, mx, tickx, ticky, label_tickx, label_ticky, title, norm='', yl=''):
    col_map=plt.get_cmap(colorcode)
    col_map.set_bad(color='gray')
    fig, ax = plt.subplots(figsize=(4,4))
    if norm!='':
      ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map, norm=norm)
    else:
      ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map)
    ax.set_xticks(tickx)
    ax.set_xticklabels(label_tickx, rotation='vertical')
    ax.set_yticks(ticky)
    ax.set_yticklabels(label_ticky)
    cbar1=fig.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=mi, vmax=mx),cmap=col_map), ax=ax)
    plt.xlabel('Phi')
    plt.ylabel(yl)
    plt.title(title)
    return ax

def plot_with_scale_bis(matrix, colorcode, mi, mx, tickx, ticky, label_tickx, label_ticky, title, yl):
    col_map=plt.get_cmap(colorcode)
    col_map.set_bad(color='gray')
    fig, ax = plt.subplots(figsize=(4,4))
    ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map, vmin=mi, vmax=mx)
    ax.set_xticks(tickx)
    ax.set_xticklabels(label_tickx, rotation='vertical')
    ax.set_yticks(ticky)
    ax.set_yticklabels(label_ticky)
    cbar1=fig.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=mi, vmax=mx),cmap=col_map), ax=ax)
    plt.xlabel('Phi')
    plt.ylabel(yl)
    plt.title(title)
    return ax

# fourier transform
def fft_data(x,dt):
    x=np.array(x)
    t=np.arange(len(x))
    fft_x=np.fft.fft(a=x)
    freq = np.fft.fftfreq(t.shape[-1])
    module_fft=pow(fft_x.real,2)+pow(fft_x.imag,2)
    freq0 = freq*(1/dt)
    freq1=freq0[1:int(len(freq0)/2)]
    module_fft0=module_fft[1:int(len(freq)/2)]
    return freq1, module_fft0

# nitrogen quotas functions
def Q_diatom(V):
    Qc=pow(10,-0.541 + 0.811*mt.log10(V))
    Qc_micro=Qc*1e-6
    Qc_micromol=Qc_micro/12
    Qn_micromol=Qc_micromol*16/106
    return Qn_micromol
def Q_eukaryotes(V):
    Qc=pow(10,-0.665 + 0.939*mt.log10(V))
    Qc_micro=Qc*1e-6
    Qc_micromol=Qc_micro/12
    Qn_micromol=Qc_micromol*16/106
    return Qn_micromol
def Q_cyanobacteria(V):
    dc=470
    Qc=dc*V
    Qc_micro=Qc*1e-9
    Qc_micromol=Qc_micro/12
    Qn_micromol=Qc_micromol*16/106
    return Qn_micromol
def Q_grazer(r):
    V=4/3*3.14159*pow(r,3)
    Qz0=pow(10,-0.547+0.9*mt.log10(V))
    Qz1=Qz0*1e-6
    Qz2=Qz1/12
    Qz3=Qz2*16/106
    return(Qz3)

def Q_virus(r):
    Na=6.022140857*1e23
    Qn=(1e6/Na)*(16*(r-2.5)**3+36*(7.5*r**2-18.75*r+15.63))
    return Qn

# encounter functions
def phiz_encounter(rh, rz, Te):
    Tk=273.15+Te
    A=1.856*1e-11
    B=4209
    C=0.04527
    D=-3.376*1e-5
    visc=A*mt.exp(B/Tk+C*Tk+D*Tk*Tk)  
    visc=visc*1e-3 #viscosity Pa.s: Pa: kg.m-1.s-2
    K= 1.380649*1e-23 #Boltzmann constant J.K-1 J=kg.m2.s-2
    pi=3.14159265359
    D_prey=K*Tk/(6*pi*visc*rh*1e-6) # host diffusion
    u_host=pow(10, 0.4 + 0.8*mt.log10(2*rh*1e-4))*1e-2
    u_zoop=pow(10,0.4 + 0.8*mt.log10(2*rz*1e-4))*1e-2
    phiz=pi*pow(3*rz*1e-6+rh*1e-6,2)*pow(u_host**2+u_zoop**2, 1/2)+4*pi*D_prey*(rz*1e-6+rh*1e-6)
    return(phiz)

def phivs_encounter(virus_radius, Host_volume, Te,indice):
    Tk=273.15+Te
    A=1.856*1e-11
    B=4209
    C=0.04527
    D=-3.376*1e-5
    visc=A*mt.exp(B/Tk+C*Tk+D*Tk*Tk)
    visc=visc*1e-3 #viscosity Pa.s: Pa: kg.m-1.s-2
    pi=3.14159265359
    Host_radius=pow(3*Host_volume/(4*pi),1/3)
    Tk=273.15+Te # water temp K
    K= 1.380649*1e-23 #Boltzmann constant J.K-1 J=kg.m2.s-2
    D_pred=K*Tk/(6*pi*visc*virus_radius*1e-9) # virus radius in m
    D_prey=0 # neglecting host diffusion

    if indice in [1,2,3]:
        u_host=pow( 10,0.4 + 0.8*mt.log10(2*Host_radius*1e-4))*1e-2 # host speed in cm.s-1 then m.s-1
    else:
        u_host=1e-6

    phi_e_a=pi*u_host*(Host_radius*1e-6+virus_radius*1e-9)**2+4*pi*(D_pred+D_prey)*(Host_radius*1e-6+virus_radius*1e-9)
    phi_e=4*pi*Host_radius*1e-6*D_pred*0.5*(1+pow(1+2*u_host*Host_radius*1e-6/D_pred, 1/3))
    phi_e_a=phi_e_a*86400*1000
    phi_e=phi_e*86400*1000

    return phi_e, phi_e_a

# load a vector
def load_vector(name, sep):
    with open(name, 'r') as f:
        line =f.readline()
        line =line.rstrip()
        l = [np.float64(num) for num in line.split(sep)]
    f.close()
    return(l)

# load a matrix
def load_matrix(name, sep):
    with open(name, 'r') as f:
        l = [[float(num) for num in line.rstrip().split(sep)] for line in f]
    f.close()
    return(l)

# write a matrix
def write_matrix(mat,name, sep):
    with open(name, 'w') as f:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                f.write(str(mat[i,j]))
                f.write(sep)
            f.write('\n')
    f.close()
# write a vector
def write_vector(vec,name, sep):
    with open(name, 'w') as f:
        for i in range(len(vec)):
            f.write(str(vec[i]))
            f.write(sep)
        f.write('\n')
# find low and high peaks from a time series
def find_low_and_high(vec):
    vec_peaks_h=find_peaks(vec)
    vec_peaks_l=find_peaks(-vec)
    res=vec_peaks_h[0]
    res_l=vec_peaks_l[0]
    if len(res)>0 and len(res_l)>0:
        vec_peak=res[0]
        if len(res)>1:
            vec_peak_b=res[1]
        else:
            vec_peak_b=float("NaN")
        vec_peak_l=res_l[0]
        if vec_peak_l<vec_peak and len(res_l)>1:
            vec_peak_l=res_l[1]
    else:
        vec_peak=float("NaN")
        vec_peak_l=float("NaN")
        vec_peak_b=float("NaN")
    return vec_peak, vec_peak_l, vec_peak_b
    
# euclidean distance between 2 vectors
def euclidean_distance(v1, v2):
    ed=mt.sqrt(sum([((v1[i]-v2[i])*100/v1[i])**2 for i in range(len(v1))]))
    return ed
    
# abosulte error between two vectors
def absolute_error(v1,v2):
    abs_err=sum([abs((v1[i]-v2[i])*100/v1[i]) for i in range(len(v1))])
    return abs_err
    
# draw borders from a matricx filled with 1 and 0: draws border around the ones
def draw_borders(matrix, lwd, edge_offset, colb, ax):
    # Draw borders
    rows, cols = matrix.shape
    for i in range(rows):
        for j in range(cols):
            if matrix[i, j] == 1:
                x = j - 0.5
                y = i - 0.5

                # Top
                if i == 0 or matrix[i - 1, j] == 0:
                    y0 = y + edge_offset if i == 0 else y
                    ax.plot([x, x + 1], [y0, y0], color=colb, linewidth=lwd)
                # Bottom
                if i == rows - 1 or matrix[i + 1, j] == 0:
                    y1 = y + 1 - edge_offset if i == rows - 1 else y + 1
                    ax.plot([x, x + 1], [y1, y1], color=colb, linewidth=lwd)
                # Left
                if j == 0 or matrix[i, j - 1] == 0:
                    x0 = x + edge_offset if j == 0 else x
                    ax.plot([x0, x0], [y, y + 1], color=colb, linewidth=lwd)
                # Right
                if j == cols - 1 or matrix[i, j + 1] == 0:
                    x1 = x + 1 - edge_offset if j == cols - 1 else x + 1
                    ax.plot([x1, x1], [y, y + 1], color=colb, linewidth=lwd)


