import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from numpy import inf
from generic_functions import *


def dSVZ(St,Vt, Zt, dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC):
    dS=dt*( (mu-d-eps*phi*Vt/Qv-St/CC)*St-phiz*Zt*St)
    dV=dt*((beta*phi*St*Vt/Qp)-(m+eps*phi*St*1/Qp+m2*Vt)*Vt)#
    dZ=dt*((phiz*eps_z*(St)-dz-dz2*Zt)*Zt)
    mZ=phiz*Zt
    mV=eps*phi*Vt/Qv
    return dS, dV, dZ, mZ, mV

def rk4_SVZ(St, Vt, Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC):
    k1S,k1V, k1Z, m1Z, m1V=dSVZ(St,Vt, Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC)
    k2S, k2V, k2Z, m2Z, m2V=dSVZ(St+0.5*k1S, Vt+0.5*k1V, Zt+0.5*k1Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC)
    k3S, k3V, k3Z, m3Z, m3V=dSVZ(St+0.5*k2S,Vt+0.5*k2V, Zt+0.5*k2Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC)
    k4S, k4V, k4Z, m4Z, m4V=dSVZ(St+k3S, Vt+k3V, Zt+k3Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC)
    Stn=St+(1.0/6.0)*(k1S + 2*k2S + 2*k3S + k4S)
    Vtn=Vt+(1.0/6.0)*(k1V + 2*k2V + 2*k3V + k4V)
    Ztn=Zt+(1.0/6.0)*(k1Z + 2*k2Z + 2*k3Z + k4Z)
    mZ=(1.0/6.0)*(m1Z+2*m2Z+2*m3Z+m4Z)
    mV=(1.0/6.0)*(m1V+2*m2V+2*m3V+m4V)
    return Stn, Vtn, Ztn, mZ, mV

def simulation_SVZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC, dt, ndays, init_conditions):
    #inittial conditions
    S0=init_conditions[0]
    V0=init_conditions[1]
    Z0=init_conditions[2]

    S=[S0]
    V=[V0]
    Z=[Z0]

    mVs=[0]
    mZs=[0]

    qm=0

    time_sim=range(round(ndays/dt))


    for t in time_sim:
        #print(t)
        St=S[t]
        Vt=V[t]
        Zt=Z[t]

        Stn, Vtn, Ztn, mZ, mV=rk4_SVZ(St,  Vt, Zt, dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC)
        Stn=max([0,Stn])
        Vtn=max([0,Vtn])
        Ztn=max([0,Ztn])


        S.append(Stn)
        V.append(Vtn)
        Z.append(Ztn)

        mZs.append(mZ)
        mVs.append(mV)
    return S, V, Z, mZs, mVs

def make_plots(S, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    S=S[i1:i2]
    V=V[i1:i2]
    Z=Z[i1:i2]
    time_sim=range(i1, i2)
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, I[0:len(I)-1])
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, V, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, Z, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim()
    plt.ylim((0,top))
    plt.title(tit)
    leg_vec=['Susceptible',  'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp)
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    #plt.ylim((0,top))
    plt.yscale('log')
    leg_vec=['Susceptible', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    mxs=np.array([np.max(np.array(S[i1:i2])/Qp),  np.max(np.array(V[i1:i2])/Qv), np.max(np.array(Z[i1:i2])/Qz)])
    maxi=10*np.max(mxs)
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp)
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    bottom, top = plt.ylim()
    #mini=max([bottom, 1e-3])
    mini=1e-1
    maxi=2e10
    plt.ylim((mini,maxi))
    leg_vec=['Susceptible', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp)
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    mini=1
    maxi=5e10
    plt.ylim((mini,maxi))
    leg_vec=['Susceptible', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()


    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S)/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, I[0:len(I)-1])
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, np.array(V)/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z)/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim()
    plt.yscale('log')
    leg_vec=['Susceptible',  'Virus', 'Zooplankton']
    plt.title(tit)
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S, color='#2ca02c')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim()
    plt.ylim((0,top))
    leg_vec=['Susceptible']
    plt.title(tit)
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S)/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, np.array(V[0:len(V)-1])/Qv, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    plt.ylim((0,top))
    leg_vec=['Susceptible']
    plt.title(tit)
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()


    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, Z, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim()
    plt.ylim((0,top))
    leg_vec=['Zooplankton']
    plt.title(tit)
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(V), color='#1f77b4')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmol.L-1)')
    leg_vec=[ 'Virus']
    bottom, top = plt.ylim()
    plt.ylim((0,top))
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(V)/Qv, color='#1f77b4')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    leg_vec=[ 'Virus']
    bottom, top = plt.ylim()
    #plt.ylim((0,top))
    plt.yscale('log')
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

#def fft_data(x,dt):
#    x=np.array(x)
#    t=np.arange(len(x))
#    t=t
#    fft_x=np.fft.fft(a=x)
#    freq = np.fft.fftfreq(t.shape[-1])
#    module_fft=pow(fft_x.real,2)+pow(fft_x.imag,2)
#    freq0 = freq*(1/dt)
#    freq1=freq0[1:int(len(freq0)/2)]
#    module_fft0=module_fft[1:int(len(freq)/2)]
#    return freq1, module_fft0


def plot_with_scale(matrix, colorcode, mi, mx, tickx, ticky, label_tickx, label_ticky, title, norm=''):
    col_map=plt.get_cmap(colorcode)
    col_map.set_bad(color='gray')
    fig, ax = plt.subplots(figsize=(4,2))
    if norm!='':
        ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map, norm=norm)
    else:
        ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map)
    ax.set_xticks(tickx)
    ax.set_xticklabels(label_tickx, rotation='vertical')
    ax.set_yticks(ticky)
    ax.set_yticklabels(label_ticky)
    #absmax=max([np.nanmin(np.array(final_ZV_ratio)), np.nanmax(np.array(final_ZV_ratio))])
    cbar1=fig.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=mi,vmax=mx),cmap=col_map), ax=ax)
    plt.xlabel('Phi')
    #plt.ylabel('Latent period')
    plt.title(title)

def plot_with_scale_bis(matrix, colorcode, mi, mx, tickx, ticky, label_tickx, label_ticky, title, norm=''):
    col_map=plt.get_cmap(colorcode)
    col_map.set_bad(color='gray')
    fig, ax = plt.subplots(figsize=(4,2))
    if norm!='':
        ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map, norm=norm)
    else:
        ax.imshow(np.flipud(np.transpose(matrix)), cmap=col_map, vmin=mi, vmax=mx)
    ax.set_xticks(tickx)
    ax.set_xticklabels(label_tickx, rotation='vertical')
    ax.set_yticks(ticky)
    ax.set_yticklabels(label_ticky)
    #absmax=max([np.nanmin(np.array(final_ZV_ratio)), np.nanmax(np.array(final_ZV_ratio))])
    cbar1=fig.colorbar(cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=mi,
                                                                   vmax=mx),
                                  cmap=col_map), ax=ax)
    plt.xlabel('Phi')
    #plt.ylabel('Latent period')
    plt.title(title)

# Reversing a list
#def Reverse(lst):
#    new_lst = lst[::-1]
#    return new_lst
#def Q_diatom(V):
#    Qc=pow(10,-0.541 + 0.811*mt.log10(V))
#    Qc_micro=Qc*1e-6
#    Qc_micromol=Qc_micro/12
#    Qn_micromol=Qc_micromol*16/106
#    return Qn_micromol

#def Q_eukaryotes(V):
#    Qc=pow(10,-0.665 + 0.939*mt.log10(V))
#    Qc_micro=Qc*1e-6
#    Qc_micromol=Qc_micro/12
#    Qn_micromol=Qc_micromol*16/106
#    return Qn_micromol
#def Q_cyanobacteria(V):
#    dc=470
#    Qc=dc*V
#    Qc_micro=Qc*1e-9
#    Qc_micromol=Qc_micro/12
#    Qn_micromol=Qc_micromol*16/106
#    return Qn_micromol
#def Q_grazer(r):
#    V=4/3*3.14159*pow(r,3)
#    Qz0=pow(10,-0.547+0.9*mt.log10(V))
#    Qz1=Qz0*1e-6
#    Qz2=Qz1/12
#    Qz3=Qz2*16/106
#    return(Qz3)

#def Q_virus(r):
#    Na=6.022140857*1e23
#    Qn=(1e6/Na)*(16*(r-2.5)**3+36*(7.5*r**2-18.75*r+15.63))
#    return Qn

#def load_vector(name, sep):
#    with open(name, 'r') as f:
#        line =f.readline()
#        line =line.rstrip()
#        l = [np.float64(num) for num in line.split(sep)]
#    f.close()
#    return(l)

#def write_vector(vec,name, sep):
#    with open(name, 'w') as f:
#        for i in range(len(vec)):
#            f.write(str(vec[i]))
#            f.write(sep)
#        f.write('\n')
#    f.close()

#def write_matrix(mat,name, sep):
#    with open(name, 'w') as f:
#        for i in range(mat.shape[0]):
#            for j in range(mat.shape[1]):
#                f.write(str(mat[i,j]))
#                f.write(sep)
#            f.write('\n')
#    f.close()

