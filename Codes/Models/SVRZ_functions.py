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
def dSVRZ(St,Vt,Rt, Zt, dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur):
    dS=dt*( (mu-d-eps*phi*Vt/Qv-(St+Rt)/CC)*St-phiz*Zt*St-mu*eps_r*St+eps_lr*mur*Rt)
    dV=dt*((beta*phi*eps*St*Vt/Qp)+(beta*phir*epsr*Rt*Vt/Qp)-(m+phi*St*1/Qp+phir*Rt/Qp+m2*Vt)*Vt)#
    dR=dt*( (mur-d-epsr*phir*Vt/Qv-(St+Rt)/CC)*Rt-phiz*Zt*Rt-mur*eps_lr*Rt+eps_r*mu*St)
    dZ=dt*((phiz*eps_z*(St+Rt)-dz-dz2*Zt)*Zt)
    mZ=phiz*Zt
    mV=eps*phi*Vt/Qv*(St/(St+Rt))+eps*phir*Vt/Qv*(Rt/(St+Rt))
    return dS, dV, dR,dZ,mZ,mV

# runge kutta 4
def rk4_SVRZ(St, Vt, Rt,Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur):
    k1S, k1V, k1R,k1Z,  m1Z, m1V=dSVRZ(St,Vt,Rt, Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur)
    k2S, k2V, k2R,k2Z, m2Z, m2V=dSVRZ(St+0.5*k1S, Vt+0.5*k1V,Rt+0.5*k1R, Zt+0.5*k1Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur)
    k3S, k3V, k3R,k3Z, m3Z, m3V=dSVRZ(St+0.5*k2S,Vt+0.5*k2V,Rt+0.5*k2R, Zt+0.5*k2Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur)
    k4S, k4V, k4R,k4Z, m4Z, m4V=dSVRZ(St+k3S, Vt+k3V, Rt+k3R,Zt+k3Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur)
    Stn=St+(1.0/6.0)*(k1S + 2*k2S + 2*k3S + k4S)
    Vtn=Vt+(1.0/6.0)*(k1V + 2*k2V + 2*k3V + k4V)
    Rtn=Rt+(1.0/6.0)*(k1R + 2*k2R + 2*k3R + k4R)
    Ztn=Zt+(1.0/6.0)*(k1Z + 2*k2Z + 2*k3Z + k4Z)

    mZ=(1.0/6.0)*(m1Z + 2*m2Z + 2*m3Z + m4Z)
    mV=(1.0/6.0)*(m1V + 2*m2V + 2*m3V + m4V)
    return Stn,Vtn,Rtn,Ztn, mZ, mV 

# simulation
def simulation_SVRZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur, dt, ndays, init_conditions):
    #inittial conditions
    S0=init_conditions[0]
    V0=init_conditions[1]
    R0=init_conditions[2]
    Z0=init_conditions[3]

    S=[S0]
    V=[V0]
    R=[R0]
    Z=[Z0]

    mZs=[0]
    mVs=[0]

    qm=0

    time_sim=range(round(ndays/dt))


    for t in time_sim:
        St=S[t]
        Vt=V[t]
        Rt=R[t]
        Zt=Z[t]

        Stn, Vtn,Rtn, Ztn, mZ, mV=rk4_SVRZ(St,  Vt,Rt, Zt, dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC, eps_r,eps_lr,phir,mur)
        Stn=max([0,Stn])
        Vtn=max([0,Vtn])
        Rtn=max([0,Rtn])
        Ztn=max([0,Ztn])


        S.append(Stn)
        V.append(Vtn)
        R.append(Rtn)
        Z.append(Ztn)

        mZs.append(mZ)
        mVs.append(mV)
    return S, V, R, Z, mZs, mVs

# plot time series 1
def make_1_plot(S, V, R,Z,i1,i2, title,dt, pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    time_sim=range(i1, i2)
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color='#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    plt.yscale('log')
    leg_vec=['Susceptible','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

# plot time series 2
def make_plots(S, V, R,Z, i1, i2, tit,dt,pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    S=S[i1:i2]
    V=V[i1:i2]
    R=R[i1:i2]
    Z=Z[i1:i2]
    time_sim=range(i1, i2)
  # molar time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, R, color='#9467bd')
    plt.plot(np.array(time_sim)*dt, V, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, Z, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    plt.ylim((1e-7,0.5))
    plt.yscale('log')
    plt.title(tit)
    leg_vec=['Susceptible', 'Resistant', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

  # count time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S)/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(R)/Qp, color='#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V)/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z)/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    mini=1
    maxi=5e10
    plt.ylim((mini,maxi))
    leg_vec=['Susceptible', 'Resistant' ,'Virus', 'Zooplankton']
    plt.title(tit)
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
