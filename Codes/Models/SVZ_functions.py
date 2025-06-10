import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from numpy import inf
from generic_functions import *

# integration
def dSVZ(St,Vt, Zt, dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,
                      epso,phiz,eps_z, dz,dz2,CC):
    dS=dt*( (mu-d-eps*phi*Vt/Qv-St/CC)*St-phiz*Zt*St)
    dV=dt*((beta*phi*St*Vt/Qp)-(m+eps*phi*St*1/Qp+m2*Vt)*Vt)#
    dZ=dt*((phiz*eps_z*(St)-dz-dz2*Zt)*Zt)
    mZ=phiz*Zt
    mV=eps*phi*Vt/Qv
    return dS, dV, dZ, mZ, mV

# runge kutta 4
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

# simulation
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

# time series plot 1
def make_plots(S, V, Z, i1, i2, tit,dt,pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    S=S[i1:i2]
    V=V[i1:i2]
    Z=Z[i1:i2]
    time_sim=range(i1, i2)
    # molar time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, V, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, Z, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    plt.ylim((1e-7,0.5))
    plt.yscale('log')
    plt.title(tit)
    leg_vec=['Susceptible',  'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    # count time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
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
