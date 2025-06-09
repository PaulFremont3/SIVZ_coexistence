import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib
import math as mt
from matplotlib.backends.backend_pdf import PdfPages
import sys
from numpy import inf
import copy
import random as random
from scipy.signal import find_peaks
from generic_functions import *

#def simulation_SIVRZ_hill(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,
#                      epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,Pc,Dc, dt, ndays):
    #inittial conditions
#    S0=0.001
#    V0=0.001
#    I0=0.001
#    R0=0.001
#    Z0=0.001

#    S=[S0]
#    I=[I0]
#    V=[V0]
#    R=[R0]
#    Z=[Z0]

#    qm=0
#
#    time_sim=range(round(ndays/dt))
#

#    for t in time_sim:
#        #print(t)
#        St=S[t]
#        It=I[t]
#        Vt=V[t]
#        Rt=R[t]
#        Zt=Z[t]

        
#        l2=pow(St+Rt+It,2)
#        ldc2=pow(Dc,2)
#        lpc2=pow(Pc,2)
#        hill_V=l2/(ldc2+l2)
#        hill_Z=l2/(lpc2+l2)
        
#        ch = ( (mu-ds-eps*phi*hill_V*Vt/Qv-(St+Rt)/CC)*St +epso*It-qm*St*St+eps_lr*mur*Rt-eps_r*mu*St-phiz*hill_Z*Zt*St)*dt
        
        #print(eps*phi*Vt/Qv*St)
#        if ch<-St:
#            Stn=St-St
#        else:
#            Stn=St+ch
#        if eps*phi*hill_V*Vt*St/Qv*dt>St:
#            #print('ok')
#            Itn =  It+(mui-ds-(1/eta))*It*dt+(St-ds*St*dt)-epso*It*dt-phiz*hill_Z*Zt*It*dt
#        else:
#            Itn =  It+(mui-ds-(1/eta))*It*dt+eps*phi*hill_V*Vt*St/Qv*dt+eps*phir*hill_V*Vt*Rt/Qv*dt-epso*It*dt-phiz*hill_Z*Zt*It*dt
#
#        chv=eps*phi*hill_V*St*Vt/Qp*dt
#        if chv>Vt:
#            Vtn= Vt+( (1/eta)*beta*It*(Qv/Qp)-m*Vt)*dt-Vt
#            chvr=eps*phir*hill_V*Rt*Vt/Qp*dt
#            if chvr>Rt:
#                Vtn= Vt-Rt
#            else:
#                Vtn=Vt-chvr
#        else:
#            Vtn= Vt+( (1/eta)*beta*It*(Qv/Qp)-(m+m2*Vt+eps*phi*hill_V*St*1/Qp+eps*phir*hill_V*Rt*1/Qp)*Vt )*dt
#            
#        Rtn = Rt +  ( (mur-ds-eps*phir*hill_V*Vt/Qv-(St+Rt)/CC)*Rt+eps_r*mu*St-eps_lr*mur*Rt-qm*Rt*Rt-phiz*hill_Z*Zt*Rt)*dt
#        
#        Ztn = Zt + (phiz*eps_z*hill_Z*(St+It+Rt)-dz-dz2*Zt)*Zt*dt
#
#        Stn=max([0,Stn])
#        Itn=max([0,Itn])
#        Vtn=max([0,Vtn])
#        Rtn=max([0,Rtn])
#        Ztn=max([0,Ztn])

        #if Stn/Qp < 1:
        #    Stn=0
        #if Itn/Qp < 1:
        #    Itn=0
        #if Rtn/Qp < 1:
        #    Rtn=0
        #if Vtn/Qv < 1:
        #    Vtn=0

#        S.append(Stn)
#        I.append(Itn)
#        V.append(Vtn)
#        R.append(Rtn)
#        Z.append(Ztn)
#    return S, I, V, R, Z

def dSIVRZ(St, It,Vt,Rt, Zt,dt, mu, mui, eta, beta, phi,ds,m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha):
    dS=( (mu-ds-eps*phi*Vt/Qv-(St+Rt)/CC-phiz*Zt-eps_r*mu)*St+eps_lr*mur*Rt )*dt
    dI=( (eps*phi*Vt*St/Qv+epsr*phir*Vt*Rt/Qv)-(ds+(1/eta)+alpha*phiz*Zt)*It )*dt
    dV=( (beta/eta)*It*(Qv/Qp)-(m+phi*St/Qp+phir*Rt/Qp+m2*Vt)*Vt )*dt #
    dR=( (mur-ds-epsr*phir*Vt/Qv-(St+Rt)/CC-phiz*Zt-eps_lr*mur)*Rt+eps_r*mu*St)*dt
    dZ=(phiz*eps_z*(St+alpha*It+Rt)-dz-dz2*Zt)*Zt*dt
    mZ=phiz*Zt
    if St+Rt+It >0:
        pi=It/(St+Rt+It)
    else:
        pi=0
    mV=pi/eta #weighted mortality linked to virus
    npp=(mu-ds)*St+(mur-ds)*Rt
    return dS, dI, dV, dR, dZ, mZ, mV, npp

def rk4_SIVRZ(St, It,Vt,Rt, Zt,dt, mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha):
    k1S, k1I, k1V, k1R, k1Z, m1Z, m1V, npp1=dSIVRZ(St, It,Vt, Rt, Zt,dt,mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha)
    k2S, k2I, k2V, k2R, k2Z, m2Z, m2V, npp2=dSIVRZ(St+0.5*k1S, It+0.5*k1I,Vt+0.5*k1V,Rt+0.5*k1R, Zt+0.5*k1Z,dt,mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha)
    k3S, k3I, k3V, k3R, k3Z, m3Z, m3V, npp3=dSIVRZ(St+0.5*k2S, It+0.5*k2I,Vt+0.5*k2V, Rt+0.5*k2R, Zt+0.5*k2Z,dt,mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha)
    k4S, k4I, k4V, k4R, k4Z, m4Z, m4V, npp4=dSIVRZ(St+k3S, It+k3I, Vt+k3V, Rt+k3R,Zt+k3Z,dt,mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha)
    Stn=St+(1.0/6.0)*(k1S + 2*k2S + 2*k3S + k4S)
    Itn=It+(1.0/6.0)*(k1I + 2*k2I + 2*k3I + k4I)
    Vtn=Vt+(1.0/6.0)*(k1V + 2*k2V + 2*k3V + k4V)
    Rtn=Rt+(1.0/6.0)*(k1R + 2*k2R + 2*k3R + k4R)
    Ztn=Zt+(1.0/6.0)*(k1Z + 2*k2Z + 2*k3Z + k4Z)

    mZ=(1.0/6.0)*(m1Z + 2*m2Z + 2*m3Z + m4Z)
    mV=(1.0/6.0)*(m1V + 2*m2V + 2*m3V + m4V)
    npp=(1.0/6.0)*(npp1 + 2*npp2 + 2*npp3 + npp4)
    return Stn,Itn,Vtn,Rtn,Ztn, mZ, mV, npp

def simulation_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,dt, ndays):
    #inittial conditions
    S0=0.001
    V0=0.001
    I0=0.001
    R0=0.001
    Z0=0.001

    S=[S0]
    I=[I0]
    V=[V0]
    R=[R0]
    Z=[Z0]

    qm=0

    time_sim=range(round(ndays/dt))


    for t in time_sim:
        #print(t)
        St=S[t]
        It=I[t]
        Vt=V[t]
        Rt=R[t]
        Zt=Z[t]

        Stn = St + ( (mu-ds-phi*Vt/Qv-(St+Rt)/CC-phiz*Zt-eps_r*mu)*St+eps_lr*mur*Rt )*dt
        
        Itn = It + ( (phi*Vt*St/Qv+phir*Vt*Rt/Qv)-(ds+(1/eta)+phiz*Zt)*It )*dt
        
        Vtn = Vt + ( (beta/eta)*It*(Qv/Qp)-(m+phi*St/Qp+phir*Rt/Qp)*Vt )*dt
            
        Rtn = Rt + ( (mur-ds-phir*Vt/Qv-(St+Rt)/CC-phiz*Zt-eps_lr*mur)*Rt+eps_r*mu*St)*dt
        
        Ztn = Zt + (phiz*eps_z*(St+It+Rt)-dz-dz2*Zt)*Zt*dt

        Stn=max([0,Stn])
        Itn=max([0,Itn])
        Vtn=max([0,Vtn])
        Rtn=max([0,Rtn])
        Ztn=max([0,Ztn])
  
        #if Stn/Qp < 1:
        #    Stn=0
       # if Itn/Qp < 1:
       #     Itn=0
       # if Rtn/Qp < 1:
       #     Rtn=0
       # if Vtn/Qv < 1:
       #     Vtn=0
       # if Ztn/Qz < 1:
       #     Ztn=0

        S.append(Stn)
        I.append(Itn)
        V.append(Vtn)
        R.append(Rtn)
        Z.append(Ztn)
    return S, I, V, R, Z


def simulation_SIVRZ_rk4(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha, dt, ndays, init):
    #inittial conditions
    S0=init[0]
    V0=init[1]
    I0=init[2]
    R0=init[3]
    Z0=init[4]

    S=[S0]
    I=[I0]
    V=[V0]
    R=[R0]
    Z=[Z0]

    mZs=[0]
    mVs=[0]
    npps=[0]

    time_sim=range(round(ndays/dt))

    for t in time_sim:
        #print(t)
        St=S[t]
        It=I[t]
        Vt=V[t]
        Rt=R[t]
        Zt=Z[t]

        Stn, Itn, Vtn, Rtn, Ztn, mZ, mV, npp=rk4_SIVRZ(St, It,Vt,Rt, Zt,dt, mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,eps,epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha)

        Stn=max([0,Stn])
        Itn=max([0,Itn])
        Vtn=max([0,Vtn])
        Rtn=max([0,Rtn])
        Ztn=max([0,Ztn])

        S.append(Stn)
        I.append(Itn)
        V.append(Vtn)
        R.append(Rtn)
        Z.append(Ztn)

        mZs.append(mZ)
        mVs.append(mV)
        npps.append(npp)
    return S, I, V, R, Z, mZs, mVs, npps

#def equilibrium1_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
#    lp=eta
#    gr=phiz
#
#    A=beta*Qv/(lp*m2*Qp)
#    a=mur-mu
#    b=mt.sqrt(A)*(epsr*phir-eps*phi)/Qv
#    try:
#        I_star_sqrt=a/b
#    except:
#        I_star_sqrt=0
#    if I_star_sqrt>0:
#        I_star=I_star_sqrt**2
#    try:
#        c=mur-ds-epsr*phir*mt.sqrt(A)*a/b+gr*dz/dz2-gr*gr*eps_z*I_star/dz2
#    except ZeroDivisionError:
#        c=0
#    d=1/CC+gr*gr*eps_z/dz2

#    try:
#        A=1/eta+ds-epsr*phir*mt.sqrt(A)*(c/d)+gr*gr*eps_z/dz2*(a/b)*(alpha)-gr*dz/dz2
        #A=1/eta+ds-beta*phir/(m*eta*Qp)*(c/d)
#    except ZeroDivisionError:
#        A=0
#    B=-b-gr*gr*eps_z/dz2
#    if B!=0:
#       S_star=A/B
#       R_star=c/d-S_star
#       RS_star=c/d
#       V_star=mt.sqrt((m+phi*S_star/Qp+phir*R_star/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star)/(2*m2)-m/(2*m2)-phi*S_star/(Qp*2*m2)-phir*R_star/(Qp*2*m2)
#       Z_star=(eps_z*gr*(I_star+S_star+R_star)-dz)/dz2
#    else:
#       S_star=0
#       R_star=0
#       V_star=0
#       Z_star=0
#    return S_star, I_star, V_star,R_star, Z_star, RS_star

def equilibrium_SVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr,
                      epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur):
    a=beta*eps*phi/(Qp*m2)-phi/(Qp*m2)
    b=beta*epsr*phir/(Qp*m2)-phir/(Qp*m2)
    c=-m/m2

    A=-eps*phi/Qv*a-1/CC-(phiz**2)*eps_z/dz2
    B=-eps*phi/Qv*b-1/CC-(phiz**2)*eps_z/dz2
    C=mu-ds-eps*phi/Qv*c+phiz*dz/dz2

    A1=-epsr*phir/Qv*a-1/CC-(phiz**2)*eps_z/dz2
    B1=-epsr*phir/Qv*b-1/CC-(phiz**2)*eps_z/dz2
    C1=mur-ds-epsr*phir/Qv*c+phiz*dz/dz2

    R_star=-(C-C1*A/A1)/(B-B1*A/A1)
    S_star=(-B1*R_star-C1)/A1

    V_star=a*S_star+b*R_star+c
    Z_star=(eps_z*phiz*(S_star+R_star)-dz)/dz2

    return S_star, R_star, V_star, Z_star

def equilibrium1_SIVRZ_m2(mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
    g=phiz

    # first approximation SVRZ
    S_star_SVRZ, R_star_SVRZ, V_star_SVRZ, Z_star_SVRZ=equilibrium_SVRZ(mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr,epso,phiz,eps_z, dz,dz2,CC,eps_r,eps_lr,phir,mur)
    
    a1=-g*g*eps_z/dz2
    b1=-ds-1/lp-(g**2)*eps_z*(S_star_SVRZ+R_star_SVRZ)/dz2+g*dz/dz2
    c1=eps*phi*S_star_SVRZ*V_star_SVRZ/Qv+epsr*phir*R_star_SVRZ*V_star_SVRZ/Qv

    try:
        if a1!=0:
            I_star0=(-b1-mt.sqrt(b1**2-4*a1*c1))/(2*a1)
            I_star1=copy.deepcopy(I_star0)
        else:
            I_star0=-c1/b1
    except:
        I_star0=0

    if I_star0>0:
        # continue refining the equilibrium
        S_star0=S_star_SVRZ
        R_star0=R_star_SVRZ
        for i in range(10):
            try:
                V_star0= (mur-mu)/((epsr*phir-eps*phi)/Qv) #mt.sqrt((m+phi*S_star0/Qp+phir*R_star0/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star0)/(2*m2)-m/(2*m2)-phi*S_star0/(Qp*2*m2)-phir*R_star0/(Qp*2*m2)
            except:
                V_star0=0
                I_star0=0
                break
            a1=-g*g*eps_z/dz2
            b1=-ds-1/lp-(g**2)*eps_z*(S_star0+R_star0)/dz2+g*dz/dz2
            c1=eps*phi*S_star0*V_star0/Qv+epsr*phir*R_star0*V_star0/Qv
            try:
                if a1!=0:
                    I_star0=(-b1-mt.sqrt(b1**2-4*a1*c1))/(2*a1)
                else:
                    I_star0=-c1/b1
            except:
                I_star0=0
            #try:
            #    C=m*V_star0+m2*V_star0**2+phi/Qp*S_star0*V_star0+phir/Qp*R_star0*V_star0
            #    D=beta*Qv/(lp*Qp)
            #    I_star0=C/D
            #except:
            #    I_star0=0
            if I_star0>0:
            # exclude the point if I_star is found to be too low
                if mt.log10(I_star0)-mt.log10(Qp)>0:
                    I_star0=I_star0
                else:
                    I_star0=0
                    break
            try:
                c=mur-ds-epsr*phir/Qp*V_star0+g*dz/dz2-g*g*eps_z*I_star0**2/dz2
            except ZeroDivisionError:
                c=0
            d=1/CC+g*g*eps_z/dz2

            RS_star0=c/d
            if RS_star0>0:
                A=I_star0/lp+ds*I_star0+(g**2)*eps_z/dz2*(I_star0**2+c/d*I_star0)-g*dz/dz2*I_star0
                B=V_star0*eps*phi/Qv
                RS_star_rat0=A/B
                R_star0=(RS_star_rat0-RS_star0)/(epsr*phir/(eps*phi)-1)
                S_star0=RS_star0-R_star0

                I_star0=(m*V_star0+m2*V_star0**2+phi*S_star0*V_star0/Qp+phir*R_star0*V_star0/Qp)/(beta*Qv/(lp*Qp))

                c=mur-ds-epsr*phir/Qp*V_star0+g*dz/dz2-g*g*eps_z*I_star0/dz2
                RS_star0=c/d

                #if R_star0>RS_star0:
                #    S_star0=(RS_star0-RS_star_rat0/(epsr*phir/(eps*phi)))/(1-epsr*phir/(eps*phi))
                #    R_star0=c/d-S_star0

                #B=1/CC+eps_z*g*g/dz2
                #A=mur-ds-epsr*phir*V_star0/Qv-S_star0*B-I_star0*eps_z*g*g/dz2+g*dz/dz2
                #R_star0=A/B
            #A_s=mur-ds-epsr*phir*V_star0/Qv-R_star0*B-I_star0*eps_z*g*g/dz2+g*dz/dz2
            #S_star0=A_s/B
        #Z_star0=(eps_z*g*(I_star0+S_star0+R_star0)-dz)/dz2
    
    if I_star0<=0 or np.isnan(I_star0):
        S_star0=0
        Z_star0=0
        V_star0=0
        R_star0=0
    else:
        Z_star0=(eps_z*g*(I_star0+S_star0+R_star0)-dz)/dz2

    
    #A=beta*Qv/(lp*m2*Qp)
    #a=mur-mu
    #b=mt.sqrt(A)*(epsr*phir-eps*phi)/Qv
    #try:
    #    I_star_sqrt=a/b
    #except:
    #    I_star_sqrt=0
    #if I_star_sqrt>0:
    #    I_star0=I_star_sqrt**2
    #    S_star0, R_star0, Z_star0=S_star_SVRZ, R_star_SVRZ, Z_star_SVRZ
    #    V_star0=mt.sqrt((m+eps*phi*S_star0/Qp+epsr*phir*R_star0/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star0)/(2*m2)-m/(2*m2)-phi*S_star0/(Qp*2*m2)-phir*R_star0/(Qp*2*m2)
    #else:
    #    I_star0=0
    #    S_star0=0
    #    Z_star0=0
    #    V_star0=0

    v1=(m+eps*phi*S_star0/Qp+epsr*phir*R_star0/Qp)**2
    v2=4*beta*Qv*m2*I_star0/(lp*Qp)

    # decide which approximation to use
    if ((v2>5*v1) and I_star0>0) : #and I_star0>0): 
        case=1
    else:
        case=2

    if case==1:
        S_star, I_star, V_star,R_star, Z_star=equilibrium1_SIVRZ_m2_case1(mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
    elif case==2 and lp<=7:
        S_star=S_star0
        I_star=I_star0
        V_star=V_star0
        R_star=R_star0
        Z_star=Z_star0
    elif case==2 and lp>7:
        S_star=np.nan
        I_star=np.nan
        V_star=np.nan
        Z_star=np.nan
        R_star=np.nan
    
    return S_star, I_star, V_star,R_star, Z_star, case

def equilibrium1_SIVRZ_m2_case1(mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
    gr=phiz

    A=beta*Qv/(lp*m2*Qp)
    a=mur-mu
    b=mt.sqrt(A)*(epsr*phir-eps*phi)/Qv
    try:
        I_star_sqrt=a/b
    except:
        I_star_sqrt=0
    if I_star_sqrt>0:
        I_star=I_star_sqrt**2
    else:
        I_star=0
    try:
        c=mur-ds-mt.sqrt(A)*epsr*phir/Qp*(a/b)+gr*dz/dz2-gr*gr*eps_z*(a/b)**2/dz2
    except ZeroDivisionError:
        c=0
    d=1/CC+gr*gr*eps_z/dz2
    
    RS_star=c/d
    V_star=mt.sqrt(A)*mt.sqrt(I_star)
    if RS_star>0:
        for i in range(2):
            A=I_star/lp+ds*I_star+(gr**2)*eps_z/dz2*(I_star**2+c/d*I_star)-gr*dz/dz2*I_star
            B=V_star*eps*phi/Qv
            RS_star_rat=A/B

            R_star=(RS_star_rat-RS_star)/(epsr*phir/(eps*phi)-1)
            S_star=RS_star-R_star
            
            C=m*V_star+m2*V_star**2+phi*S_star*V_star/Qp+phir*R_star*V_star/Qp
            D=beta*Qv/(lp*Qp)
            I_star=(C/D)
           
            try:
                c=mur-ds-mt.sqrt(A)*epsr*phir/Qp*(a/b)+gr*dz/dz2-gr*gr*eps_z*I_star/dz2
            except (ZeroDivisionError, ValueError):
                c=0
            RS_star=c/d
        Z_star=(eps_z*gr*(I_star+S_star+R_star)-dz)/dz2
    else:
        I_star=0
        S_star=0
        R_star=0
        V_star=0
        Z_star=0
    return S_star, I_star, V_star,R_star, Z_star

def equilibrium2_SIVRZ_m2(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alph):
    g=phiz

    # decide which approximation
    g=phiz
    a0=mu-d+eps_z*g*(dz/dz2)+beta*eps*phi/Qp*(m/m2)
    b0=1/CC+g*g*eps_z/dz2+beta*phi*phi*eps*eps/(Qp*m2*Qv)
    S_star_SVZ=a0/b0
    V_star_SVZ=(beta*phi/Qp*S_star_SVZ-m)/m2

    a1=-g*g*eps_z/dz2
    b1=-d-1/lp-(g**2)*eps_z*S_star_SVZ/dz2+g*dz/dz2
    c1=eps*phi*S_star_SVZ*V_star_SVZ/Qv

    try:
        if a1!=0:
            I_star0=(-b1-mt.sqrt(b1**2-4*a1*c1))/(2*a1)
            I_star1=copy.deepcopy(I_star0)
        else:
            I_star0=-c1/b1
    except:
        I_star0=0

    if I_star0>0:
        # continue refining the equilibrium
        for i in range(1000):
            S_star0=S_star_SVZ
            Z_star0=(eps_z*g*(I_star0+S_star0)-dz)/dz2
            if S_star0>0:
                V_star0= mt.sqrt((m+phi*S_star0/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star0)/(2*m2)-m/(2*m2)-phi*S_star0/(Qp*2*m2)
            else:
                V_star0=0

            a1=-g*g*eps_z/dz2
            b1=-d-1/lp-(g**2)*eps_z*S_star0/dz2+g*dz/dz2
            c1=eps*phi*S_star0*V_star0/Qv
            try:
                if a1!=0:
                    I_star0=(-b1-mt.sqrt(b1**2-4*a1*c1))/(2*a1)
                else:
                    I_star0=-c1/b1
            except:
                I_star0=0
            if I_star0>0:
                if mt.log10(I_star0)-mt.log10(Qp)>0:
                    I_star0=I_star0
                else:
                    I_star0=0
                    break
    else:
        S_star0=0
        Z_star0=0
        V_star0=0

    v1=(m+phi*S_star_SVZ/Qp)**2
    v2=4*beta*Qv*m2*I_star0/(lp*Qp)
    rat=v1/v2

    if ((v2>10*v1) and I_star0>0): 
        case=1
    else:
        case=2
    
    if case==1 and I_star0>0:
        A=beta*Qv/(lp*Qp*m2)

        denom=1/CC+g*g*eps_z/dz2-eps*phi*phi/(Qv*Qp*2*m2)

        B=(eps*phi/Qv*mt.sqrt(A))/denom
        C=(mu-d+g*dz/dz2)/denom
        D=((g**2)*eps_z/dz2)/denom
        E=-d-1/lp+g*dz/dz2
        F=eps*phi*phi/(Qv*Qp*2*m2)

        
        #a=-eps_z*g*g/dz2*(1-D)+(D**2)*eps*F
        #b=eps_z*g*g*B/dz2-D*phi/Qv*mt.sqrt(A)+2*D*B*eps*F
        #c=-eps*phi/Qv*B*mt.sqrt(A)+E-C*eps_z*g**2/dz2+eps*F-2*D*C*eps*F
        #d1=eps*phi/Qv*C*mt.sqrt(A)-eps*F*2*C*B
        #e=(eps*phi**2)*(C**2)/(Qv*Qp*2*m2)
        
        a=-eps_z*g*g/dz2*(1-D)+(D**2)*phi*phi/(Qv*Qp*2*m2)
        b=eps_z*g*g*B/dz2-D*phi/Qv*mt.sqrt(A)-2*D*B*F
        c=-eps*phi/Qv*B*mt.sqrt(A)+E-C*eps_z*g**2/dz2-F*(B**2)+2*D*C*F+m/(2*m2)*D*(eps*phi/Qv)
        d1=eps*phi/Qv*C*mt.sqrt(A)+F*2*C*B+m/(2*m2)*B*(eps*phi/Qv)
        e=-(eps*phi**2)*(C**2)/(Qv*Qp*2*m2)-m*C/(2*m2)*(eps*phi/Qv)

        try:
            roots=np.roots([a,b,c,d1,e])
            I_star_sqrt=(-d1-mt.sqrt(d1**2-4*c*e))/(2*c)
            I_star_sqrt=0
            for r in roots:
                if r.real>0 and r.real<CC*mu:
                    I_star_sqrt=r.real
                    break
        except:
            I_star_sqrt=0
        if I_star_sqrt>=0:
            I_star=pow(I_star_sqrt,2)
            S_star=C-B*mt.sqrt(I_star)-D*I_star
            V_star=mt.sqrt((m+phi*S_star/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star)/(2*m2)-m/(2*m2)-phi*S_star/(Qp*2*m2) #mt.sqrt(beta*Qv/(lp*m2*Qp))*mt.sqrt(I_star)-m/(2*m2)-phi*S_star/(Qp*2*m2)
            Z_star=(eps_z*g*(I_star+S_star)-dz)/dz2
        else:
            I_star=0
            S_star=0
            V_star=0
            Z_star=0
    elif case==1 and I_star0<0:
        I_star=0
        S_star=0
        V_star=0
        Z_star=0
    elif case==2 and lp<=10:
        S_star=S_star0
        I_star=I_star0
        V_star=V_star0
        Z_star=Z_star0
    elif case==2 and lp>10:
        S_star=np.nan
        I_star=np.nan
        V_star=np.nan
        Z_star=np.nan
    if S_star>0 and V_star>0:
        R_star=-mu*eps_r*S_star/(mur-d-S_star/CC-g*Z_star-epsr*phir*V_star/Qv-mur*eps_r)
    else:
        R_star=0
    return S_star, I_star, V_star,R_star, Z_star


def equilibrium1_SIVRZ(mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
    g=phiz
    # 1st eq S, I, V, R, Z
    st=np.nan
    a=mur-mu
    b=(beta/(lp*m*Qp))*(epsr*phir-eps*phi)
    try:
        I_star=a/b
    except ZeroDivisionError:
        I_star=0
    try:
        c=mur-ds-epsr*phir*beta/(lp*m*Qp)*I_star+g*dz/dz2-g*g*eps_z*alpha*I_star/dz2
    except ZeroDivisionError:
        c=0
    d=1/CC+g*g*eps_z/dz2
    
    RS_star=c/d
    V_star=(mur-mu)/((epsr*phir-eps*phi)/Qv)
    if RS_star>0:
        for i in range(2):
            A=I_star/lp+ds*I_star+(g**2)*eps_z/dz2*(I_star**2+c/d*I_star)-g*dz/dz2*I_star
            B=V_star*phi*eps/Qv
            RS_star_rat=A/B

            R_star=(RS_star_rat-RS_star)/(epsr*phir/(eps*phi)-1)
            S_star=RS_star-R_star 
            Z_star=(eps_z*g*(I_star+S_star+R_star)-dz)/dz2

            C=m*V_star+phi*S_star*V_star/Qp+phir*R_star*V_star/Qp
            D=beta*Qv/(lp*Qp)
            I_star=(C/D)

            c=mur-ds-epsr*phir/Qp*V_star+g*dz/dz2-g*g*eps_z*I_star/dz2
            RS_star=c/d

    else:
        S_star=0
        R_star=0
        V_star=0
        Z_star=0
    
    return S_star, I_star, V_star, R_star, Z_star

#def equilibrium1(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
#    gr=phiz
#    alpha=1
#    # 1st eq S, I, V, R, Z
#    st=np.nan
#    alpha=1
#    a=mur-mu
#    b=(beta/(eta*m*Qp))*(epsr*phir-eps*phi)
#    try:
#        I_star=a/b
#    except ZeroDivisionError:
#        I_star=0
#    try:
#        c=mur-ds-epsr*phir*beta/(eta*m*Qp)*I_star+gr*dz/dz2-gr*gr*eps_z*I_star/dz2
#    except ZeroDivisionError:
#        c=0
#    d=1/CC+gr*gr*eps_z/dz2
#
#    try:
#        A=1/eta+ds-beta*epsr*phir/(m*eta*Qp)*(c/d)+gr*gr*eps_z/dz2*(a/b)*(alpha)-gr*dz/dz2
#        #A=1/eta+ds-beta*phir/(m*eta*Qp)*(c/d)
#    except ZeroDivisionError:
#        A=0
#    B=beta*(eps*phi-epsr*phir)/(m*eta*Qp)-gr*gr*eps_z/dz2
#    if B!=0:
#       S_star=A/B
#       R_star=c/d-S_star
#       V_star=beta*Qv*I_star/(m*eta*Qp)
#       Z_star=(eps_z*gr*(I_star+S_star+R_star)-dz)/dz2
#    else:
#       S_star=0
#       R_star=0
#       V_star=0
#       Z_star=0
#    return S_star, I_star, V_star, R_star, Z_star

def equilibrium1_SIVRZ_bis(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
    gr=phiz
    a=mur-mu
    b=(beta/(eta*m*Qp))*(epsr*phir-eps*phi)
    try:
        c=mur-ds-epsr*phir*beta/(eta*m*Qp)*(a/b)
    except ZeroDivisionError:
        #print(a)
        #print(b)
        #print(phir)
        c=0
    d=1/CC+gr*gr*eps_z/dz2
    A=1/eta+d-epsr*beta*phir*c/(m*eta*Qp*d)
    B=(beta/(m*eta*Qp)*(eps*phi-epsr*phir))-gr*gr*eps_z/dz2
    if B!=0:
        S_star=A/B
    try:
        I_star=a/b
    except ZeroDivisionError:
        I_star=0
    R_star=c/d-S_star
    
    for i in range(2):
        V_star=beta*Qv*I_star/((m+phi*S_star/Qp+phir*R_star/Qp)*eta*Qp)
        C=m*V_star+phi*S_star*V_star/Qp+phir*R_star*V_star/Qp
        D=beta*Qv/(eta*Qp)
        I_star=(C/D)

    Z_star=(eps_z*gr*(S_star+R_star)-dz)/dz2
    return S_star, I_star, V_star, R_star, Z_star

#def equilibrium2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
#    gr=phiz
#    a0=mur-ds+gr*(dz/dz2)
#    b0=1/CC+gr*gr*eps_z/dz2
#    T0=(epsr*phir*beta)/(m*eta*Qp)+gr*gr*eps_z*alpha/dz2
#    try:
#        R_star0=(1/eta+ds+gr*gr*eps_z/dz2*(a0*alpha/T0-dz))/((epsr*phir*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2*(b0*alpha/T0-1))
#    except ZeroDivisionError:
#        R_star0=0
#    try:
#        I_star0=(a0-b0*R_star0)/T0
#    except ZeroDivisionError:
#        I_star0=0
    
#    V_star0=beta*Qv*I_star0/((m+phi*R_star0/Qp)*eta*Qp)
    #V_star0=beta*Qv*I_star0/(m*eta*Qp)

    
#    for i in range(1000):
#        V_star0=beta*Qv*I_star0/((m+phi*R_star0/Qp)*eta*Qp)
        #C=m*V_star0+phi*R_star0*V_star0/Qp
        #D=beta*Qv/(eta*Qp)
        #I_star0=(C/D)
        
#        A=-gr**2*eps_z/dz2
#        B=-1/eta-ds+gr*dz/dz2-gr**2*eps_z/dz2*R_star0
#        C=epsr*phir*R_star0*V_star0/Qv
#
#        try:
#            I_star0=(-B-mt.sqrt(B**2-4*A*C))/(2*A)
#            I_star0=-C/B
#        except:
#            I_star0=0
#            R_star0=0
#            break
#
#        D=mur-ds+gr*(dz/dz2)-epsr*phir*V_star0/Qv-gr**2*eps_z/dz2*I_star0
#        E=1/CC+gr**2*eps_z/dz2
#        R_star0=D/E
#
#    for i in range(1):
#        V_star0=beta*Qv*I_star0/((m+phir*R_star0/Qp)*eta*Qp)
        #C=m*V_star0+phir*R_star0*V_star0/Qp
        #D=beta*Qv/(eta*Qp)
        #I_star0=(C/D)

#    Z_star0=(eps_z*gr*(alpha*I_star0+R_star0)-dz)/dz2
#    S_star0=-mur*eps_r*R_star0/(mu-ds-R_star0/CC-gr*Z_star0-eps*phi*V_star0/Qv-mu*eps_r)
#    return S_star0, I_star0, V_star0, R_star0, Z_star0

def equilibrium2_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha):
    gr=phiz
    a1=mu-ds+gr*(dz/dz2)
    b1=1/CC+gr*gr*eps_z/dz2
    T1=(eps*phi*beta)/(m*eta*Qp)+gr*gr*eps_z*alpha/dz2
    try:
        S_star1=(1/eta+ds+gr*gr*eps_z/dz2*(a1*alpha/T1-dz))/((eps*phi*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2*(b1*alpha/T1-1))
    except ZeroDivisionError:
        S_star1=0
    I_star1=(a1-b1*S_star1)/T1
    V_star1=beta*Qv*I_star1/((m+phi*S_star1/Qp)*eta*Qp)
    #V_star1=beta*Qv*I_star1/(m*eta*Qp)
    #for i in range(1000):
    #    V_star1=beta*Qv*I_star1/((m+phi*S_star1/Qp)*eta*Qp)
        #C=m*V_star1+phi*S_star1*V_star1/Qp
        #D=beta*Qv/(eta*Qp)
        #I_star1=(C/D)
        
    #    A=-gr**2*eps_z/dz2
    #    B=-1/eta-ds+gr*dz/dz2-gr**2*eps_z/dz2*S_star1
    #    C=eps*phi*S_star1*V_star1/Qv
    #    #I_star1=-C/B
    #    if A!=0:
    #        try:
    #            I_star1=(-B-mt.sqrt(B**2-4*A*C))/(2*A)
    #        except:
    #            I_star1=0
    #            R_star1=0
    #            break
    #    else:
    #        try:
    #            I_star1=-C/B
    #        except:
    #            I_star1=0
    #            R_star1=0
    #            break

    #    D=mu-ds+gr*(dz/dz2)-eps*phi*V_star1/Qv-gr**2*eps_z/dz2*I_star1
    #    E=1/CC+gr**2*eps_z/dz2
    #    S_star1=D/E
    
    Z_star1=(eps_z*gr*(alpha*I_star1+S_star1)-dz)/dz2
    R_star1=-mu*eps_r*S_star1/(mur-ds-S_star1/CC-gr*Z_star1-epsr*phir*V_star1/Qv-mur*eps_r)
    return S_star1, I_star1, V_star1, R_star1, Z_star1

def check_stable_states_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC):
    gr=phiz
    g=phiz
    alpha=1
    # 1st eq S, I, V, R, Z
    ast=np.nan
    ost=np.nan
    fst=np.nan
    S_star, I_star, V_star, R_star, Z_star, case=equilibrium1_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star>0 and Z_star>0 and S_star>0 and R_star>0 and I_star>0:
        #print('ok_')
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, Z_star)
        ast, ost, fst=check_eigs(eigs, eig_vecs, ast, ost)

    # 2nd eq S,I,V, 0,Z
    ast1=np.nan
    ost1=np.nan
    fst1=np.nan
    S_star1, I_star1, V_star1, dR_star1, Z_star1=equilibrium2_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star1>0 and Z_star1>0 and S_star1>0 and I_star1>0 and dR_star1>0:
        #print('ok_1')
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star1, I_star1, V_star1, dR_star1, Z_star1)
        ast1, ost1, fst1=check_eigs(eigs, eig_vecs, ast1, ost1)

    # 3rd eq 0,I,V,R,Z
    ast2=np.nan
    ost2=np.nan
    fst2=np.nan                                                 
    R_star0, I_star0, V_star0, dS_star0, Z_star0=equilibrium2_SIVRZ_m2(mur, mui, eta, beta, phir, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star0>0 and Z_star0>0 and R_star0>0 and I_star0>0 and dS_star0>0:# and dS_star0>0:
        #print('ok_2')
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star0, I_star0, V_star0, R_star0, Z_star0)
        ast2, ost2, fst2=check_eigs(eigs, eig_vecs, ast2, ost2)

    # 4th eq S,I,V,R,0
    ast3=np.nan
    ost3=np.nan
    fst3=np.nan
    S_star, I_star, V_star, R_star, Z_star, case=equilibrium1_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,0,0,0,dz2, CC, alpha)
    Z_star=0
    if V_star>0  and R_star>0 and S_star>0 and I_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, Z_star)
        ast3, ost3, fst3=check_eigs(eigs, eig_vecs, ast3, ost3)

    # 5th eq S,I,V,0,0
    ast4=np.nan
    ost4=np.nan
    fst4=np.nan
    S_star, I_star, V_star, dR_star, Z_star=equilibrium2_SIVRZ_m2(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,0,0,0,dz2, CC, alpha)
    Z_star=0
    if V_star>0  and S_star>0 and I_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, dR_star, Z_star)
        ast4, ost4, fst4=check_eigs(eigs, eig_vecs, ast4, ost4)

    # 6th eq 0,I,V,R,0
    ast5=np.nan
    ost5=np.nan
    fst5=np.nan
    R_star0, I_star0, V_star0, dS_star0, Z_star0=equilibrium2_SIVRZ_m2(mur, mui, eta, beta, phir, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,0,0,0,dz2, CC, alpha)
    Z_star0=0
    if V_star0>0 and R_star0>0 and I_star0>0 and dS_star0>0:# and dS_star0>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,dS_star0, I_star0, V_star0, R_star0, Z_star0)
        ast5, ost5, fst5=check_eigs(eigs, eig_vecs, ast5, ost5)

    # 7th eq S,0,0,0,Z
    S_star=(mu-ds+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*S_star-dz)/dz2
    dR_star=-mu*eps_r*S_star/(mur-ds-S_star/CC-g*Z_star-mur*eps_r)
    ast6=np.nan
    ost6=np.nan
    fst6=np.nan
    if  Z_star>0 and S_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star,0, 0, dR_star, Z_star)
        ast6, ost6, fst6=check_eigs(eigs, eig_vecs, ast6, ost6)

    # 8th eq 0,0,0,R,Z
    R_star=(mur-ds+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*R_star-dz)/dz2
    dS_star=-mur*eps_r*R_star/(mu-ds-R_star/CC-g*Z_star-mu*eps_r)
    ast7=np.nan
    ost7=np.nan
    fst7=np.nan
    if  Z_star>0 and R_star>0 and dS_star>0:# and dS_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star,0, 0, R_star, Z_star)
        ast7, ost7, fst7=check_eigs(eigs, eig_vecs, ast7, ost7)

    # 9th eq S,0,0,0,0
    S_star=(mu-ds)*CC
    dR_star=-mu*eps_r*S_star/(mur-ds-S_star/CC-mur*eps_r)
    ast8=np.nan
    ost8=np.nan
    fst8=np.nan
    if  S_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star,0, 0, dR_star, 0)
        ast8, ost8, fst8=check_eigs(eigs, eig_vecs, ast8, ost8)

    # 10th eq 0,0,0,R,0
    R_star=(mur-ds)*CC
    dS_star=-mur*eps_r*R_star/(mu-ds-R_star/CC-mu*eps_r)
    ast9=np.nan
    ost9=np.nan
    fst9=np.nan
    if  R_star>0 and dS_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star,0, 0, R_star, 0)
        ast9, ost9, fst9=check_eigs(eigs, eig_vecs, ast9, ost9)

    # 11th eq 0,0,0,0,0
    ast10=np.nan
    ost10=np.nan
    fst10=np.nan
    eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, 0,0, 0, 0, 0)
    ast10, ost10, fst10=check_eigs(eigs, eig_vecs, ast10, ost10)

    astable=[ast, ast1,ast2,ast3,ast4,ast5,ast6,ast7, ast8, ast9, ast10]
    ostable=[ost, ost1,ost2,ost3,ost4,ost5,ost6,ost7, ost8, ost9, ost10]
    fstable=[fst, fst1,fst2,fst3,fst4,fst5,fst6,fst7, fst8, fst9, fst10]
    return astable, ostable,fstable


def check_stable_states_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC):
    
    gr=phiz
    g=phiz
    alpha=1
    # 1st eq S, I, V, R, Z
    ast=np.nan
    ost=np.nan
    fst=np.nan
    
    #a=mur-mu
    #b=(beta/(eta*m*Qp))*(epsr*phir-eps*phi)
    #try:
    #    I_star=a/b
    #except ZeroDivisionError:
    #    I_star=0
    #try:
    #    c=mur-ds-epsr*phir*beta/(eta*m*Qp)*I_star+g*dz/dz2-g*g*eps_z*I_star/dz2
    #except ZeroDivisionError:
    #    c=0
    #d=1/CC+gr*gr*eps_z/dz2

   # try:
   #     A=1/eta+ds-epsr*beta*phir/(m*eta*Qp)*(c/d)+gr*gr*eps_z/dz2*((a/b)*(alpha)+c/d)-gr*dz/dz2#1/eta+ds-beta*phir/(m*eta*Qp)*(c/d)+gr*gr*eps_z/dz2*(a/b)*(alpha)-gr*dz/dz2
        #A=1/eta+ds-beta*phir/(m*eta*Qp)*(c/d)
   # except ZeroDivisionError:
   #     A=0
   # B=beta*(eps*phi-epsr*phir)/(m*eta*Qp)-gr*gr*eps_z/dz2
   # if B!=0:
   #    S_star=A/B
   #    R_star=c/d-S_star
   #    V_star=beta*Qv*I_star/(m*eta*Qp)
   #    Z_star=(eps_z*gr*(I_star+S_star+R_star)-dz)/dz2     
   # else:
   #    S_star=0
   #    R_star=0
   #    V_star=0
   #    Z_star=       0                                 mu, mui, lp, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha
    S_star, I_star, V_star, R_star, Z_star=equilibrium1_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star>0 and Z_star>0 and S_star>0 and R_star>0 and I_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, Z_star)
        ast, ost, fst=check_eigs(eigs, eig_vecs, ast, ost)

    # 2nd eq S,I,V, 0,Z
    ast1=np.nan
    ost1=np.nan
    fst1=np.nan
    #a1=mu-ds+gr*(dz/dz2)
    #b1=1/CC+gr*gr*eps_z/dz2
    #T1=(phi*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2
    #S_star1=(1/eta+ds+gr*gr*eps_z/dz2*(a1/T1-dz))/((phi*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2*(b1/T1-1))
    #I_star1=(a1-b1*S_star1)/T1
    #V_star_ap1=beta*Qv*I_star1/(m*eta*Qp)
    #V_star1=beta*Qv*I_star1/(eta*Qp)/(m+phi*S_star1/Qp)
    #Z_star1=(eps_z*gr*(I_star1+S_star1)-dz)/dz2
    #dR_star1 = -mu*eps_r*S_star1/(mur-ds-S_star1/CC-g*Z_star1-phir*V_star1/Qv-mur*eps_r) 
    S_star1, I_star1, V_star1, dR_star1, Z_star1=equilibrium2_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star1>0 and Z_star1>0 and S_star1>0 and I_star1>0 and dR_star1>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star1, I_star1, V_star1, dR_star1, Z_star1)
        ast1, ost1, fst1=check_eigs(eigs, eig_vecs, ast1, ost1)

    # 3rd eq 0,I,V,R,Z
    ast2=np.nan
    ost2=np.nan
    fst2=np.nan
    #a0=mur-ds+gr*(dz/dz2)
    #b0=1/CC+gr*gr*eps_z/dz2
    #T0=(phir*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2
    #R_star0=(1/eta+ds+gr*gr*eps_z/dz2*(a0/T0-dz))/((phir*beta)/(m*eta*Qp)+gr*gr*eps_z/dz2*(b0/T0-1))
    #I_star0=(a0-b0*R_star0)/T0
    #V_star_ap0=beta*Qv*I_star0/(m*eta*Qp)
    #V_star0=beta*Qv*I_star0/(eta*Qp)/(m+phir*R_star0/Qp)
    #Z_star0=(eps_z*gr*(I_star0+R_star0)-dz)/dz2
    #dS_star1 = -mur*eps_r*R_star0/(mu-ds-R_star0/CC-g*Z_star0-phi*V_star0/Qv-mu*eps_r)
    R_star0, I_star0, V_star0, dS_star0, Z_star0=equilibrium2_SIVRZ(mur, mui, eta, beta, phir, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,phiz,eps_z,dz,dz2, CC, alpha)
    if V_star0>0 and Z_star0>0 and R_star0>0 and I_star0>0 and dS_star0>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, 0, I_star0, V_star0, R_star0, Z_star0)
        ast2, ost2, fst2=check_eigs(eigs, eig_vecs, ast2, ost2)

    # 4th eq S,I,V,R,0
    ast3=np.nan
    ost3=np.nan
    fst3=np.nan
    #a=mur-mu
    #b=(beta/(eta*m*Qp))*(phir-phi)
    #try:
    #    c=mur-ds-phir*beta/(eta*m*Qp)*(a/b)
    #except ZeroDivisionError:
    #    c=0
    #d=1/CC
    #try:
    #    I_star=a/b
    #except ZeroDivisionError:
    #    I_star=0
    #try:
    #   A=1/eta+ds-beta*phir/(m*eta*Qp)*(c/d)
    #except ZeroDivisionError:
    #   A=0
    #B=beta*(phi-phir)/(m*eta*Qp)
    #if B!=0:
    #    S_star=A/B
    #    R_star=c/d-S_star
    #    V_star=beta*Qv*I_star/(m*eta*Qp)
    #    Z_star=0
    S_star, I_star, V_star, R_star, Z_star=equilibrium1_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,0,0,0,dz2, CC, alpha)
    if V_star>0 and R_star>0 and S_star>0 and I_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, R_star, 0)
        ast3, ost3, fst3=check_eigs(eigs, eig_vecs, ast3, ost3)

    # 5th eq S,I,V,0,0
    ast4=np.nan
    ost4=np.nan
    fst4=np.nan
    #S_star=(ds+1/eta)/(phi*beta/(m*eta*Qp))
    #I_star=(mu-ds-S_star/CC)/(phi*beta/(m*eta*Qp))
    #V_star_ap=beta*Qv*I_star/(m*eta*Qp)
    #V_star=beta*Qv*I_star/(eta*Qp)/(m+phi*S_star/Qp)
    #dR_star=-mu*eps_r*S_star/(mur-ds-S_star/CC-phir*V_star/Qv-mur*eps_r)
    
    S_star, I_star, V_star, dR_star, Z_star=equilibrium2_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  eps,epsr, epso,eps_r,eps_lr,phir,mur,0,0,0,dz2, CC, alpha)
    if  S_star>0 and V_star>0 and I_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, dR_star, 0)
        ast4, ost4, fst4=check_eigs(eigs, eig_vecs, ast4, ost4)

    # 6th eq 0,I,V,R,0
    ast5=np.nan
    ost5=np.nan
    fst5=np.nan
    #R_star=(ds+1/eta)/(phir*beta/(m*eta*Qp))
    #I_star=(mur-ds-R_star/CC)/(phir*beta/(m*eta*Qp))
    #V_star_ap=beta*Qv*I_star/(m*eta*Qp)
    #V_star=beta*Qv*I_star/(eta*Qp)/(m+phir*R_star/Qp)
    #dS_star = -mur*eps_r*R_star/(mu-ds-R_star/CC-phi*V_star/Qv-mu*eps_r)
    R_star, I_star, V_star, dS_star, Z_star=equilibrium2_SIVRZ(mur, mui, eta, beta, phi, ds, m,m2, Qv, Qp,Qz,  epsr,eps, epso,eps_r,eps_lr,phi,mu,0,0,0,dz2, CC, alpha)
    if  R_star>0 and V_star>0 and I_star>0 and dS_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star, I_star, V_star, R_star, 0)
        ast5, ost5, fst5=check_eigs(eigs, eig_vecs, ast5, ost5)

    # 7th eq S,0,0,0,Z
    S_star=(mu-ds+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*S_star-dz)/dz2
    dR_star=-mu*eps_r*S_star/(mur-ds-S_star/CC-g*Z_star-mur*eps_r)
    ast6=np.nan
    ost6=np.nan
    fst6=np.nan
    if  Z_star>0 and S_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star,0, 0, dR_star, Z_star)
        ast6, ost6, fst6=check_eigs(eigs, eig_vecs, ast6, ost6)

    # 8th eq 0,0,0,R,Z
    R_star=(mur-ds+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*R_star-dz)/dz2
    dS_star=-mur*eps_r*R_star/(mu-ds-R_star/CC-g*Z_star-mu*eps_r)
    ast7=np.nan
    ost7=np.nan
    fst7=np.nan
    if  Z_star>0 and R_star>0 and dS_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star,0, 0, R_star, Z_star)
        ast7, ost7, fst7=check_eigs(eigs, eig_vecs, ast7, ost7)

    # 9th eq S,0,0,0,0
    S_star=(mu-ds)*CC
    dR_star=-mu*eps_r*S_star/(mur-ds-S_star/CC-mur*eps_r)
    ast8=np.nan
    ost8=np.nan
    fst8=np.nan
    if  S_star>0 and dR_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, S_star,0, 0, dR_star, 0)
        ast8, ost8, fst8=check_eigs(eigs, eig_vecs, ast8, ost8)

    # 10th eq 0,0,0,R,0
    R_star=(mur-ds)*CC
    dS_star=-mur*eps_r*R_star/(mu-ds-R_star/CC-mu*eps_r)
    ast9=np.nan
    ost9=np.nan
    fst9=np.nan
    if  R_star>0 and dS_star>0:
        eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, dS_star,0, 0, R_star, 0)
        ast9, ost9, fst9=check_eigs(eigs, eig_vecs, ast9, ost9)
    
    # 11th eq 0,0,0,0,0
    ast10=np.nan
    ost10=np.nan
    fst10=np.nan
    eigs, eig_vecs=Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps, epsr,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, 0,0, 0, 0, 0)
    ast10, ost10, fst10=check_eigs(eigs, eig_vecs, ast10, ost10)

    astable=[ast, ast1,ast2,ast3,ast4,ast5,ast6,ast7, ast8, ast9, ast10]
    ostable=[ost, ost1,ost2,ost3,ost4,ost5,ost6,ost7, ost8, ost9, ost10]
    fstable=[fst, fst1,fst2,fst3,fst4,fst5,fst6,fst7, fst8, fst9, fst10]
    return astable, ostable,fstable

def make_1_plot_SIVRZ(S, I, V, R,Z,i1,i2, title,dt, pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    time_sim=range(i1, i2)
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    #plt.ylim((0,top))
    plt.yscale('log')
    leg_vec=['Susceptible', 'Infected','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()


def make_plots_SIVRZ(S, I, V, R, Z, i1, i2, tit,dt, pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    time_sim=range(i1,i2)
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S[i1:i2], color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, I[i1:i2], color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, R[i1:i2], color= '#9467bd')
    plt.plot(np.array(time_sim)*dt, V[i1:i2], color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, Z[i1:i2], color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim() 
    #top=min([top, 1e20])
    #plt.ylim((0,top))
    plt.ylim((1e-7,0.5))
    plt.yscale('log')
    leg_vec=['Susceptible', 'Infected','Resistant', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    bottom, top = plt.ylim()
    top=min([top, 1e20])
    plt.ylim((bottom,top))
    leg_vec=['Susceptible', 'Infected','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    #plt.show()
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    mxs=np.array([np.max(np.array(S[i1:i2])/Qp),np.max(np.array(R[i1:i2])/Qp) ,np.max(np.array(I[i1:i2])/Qp), np.max(np.array(V[i1:i2])/Qv), np.max(np.array(Z[i1:i2])/Qz)])
    maxi=10*np.max(mxs) 
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    bottom, top = plt.ylim()
    mini=1e-1#max([bottom, 1e-3])
    top=2e10#min([maxi, 1e20])
    plt.ylim((mini,top))
    leg_vec=['Susceptible', 'Infected','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    #plt.show()
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    mini=1
    maxi=5e10
    plt.ylim((mini,maxi))
    leg_vec=['Susceptible', 'Infected','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    #plt.show()
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    #bottom, top = plt.ylim()
    #plt.ylim((0,top))
    plt.yscale('log')
    bottom, top = plt.ylim()
    top=min([top, 1e20])
    plt.ylim((1,top))
    leg_vec=['Susceptible', 'Infected','Resistant' ,'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    #plt.show()
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)
    
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S[i1:i2], color='#2ca02c')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim()
    top=min([top, 1e20])
    plt.ylim((0,top))
    leg_vec=['Susceptible']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    fig, ax = plt.subplots(figsize=(4,4))    
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    #plt.plot(np.array(time_sim)*dt, np.array(V[0:len(V)-1])/Qv, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    top=min([top, 1e20])
    plt.ylim((0,top))
    leg_vec=['Susceptible']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(R[i1:i2])/Qp, color= '#9467bd')
    #plt.plot(np.array(time_sim)*dt, np.array(V[0:len(V)-1])/Qv, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    top=min([top, 1e20])
    plt.ylim((0,top))
    leg_vec=['Susceptible', 'Infected', 'Resistant']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()
    #print(top)

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, Z[i1:i2], color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (mmolN.L-1)')
    bottom, top = plt.ylim() 
    top=min([top, 1e20])
    plt.ylim((0,top))
    #print(top)
    leg_vec=['Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    #print(top)
    top=min([top, 1e20])
    #plt.ylim((0,top))
    plt.yscale('log')
    leg_vec=['Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    leg_vec=[ 'Virus']
    top=min([top, 1e20])
    bottom, top = plt.ylim() 
    #print(top)
    plt.ylim((0,top))
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(tit)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()


def Jacobian_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epsr, epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC, Ss, Is, Vs, Rs, Zs):
    Jac=np.array([[mu-ds-eps*phi*Vs/(Qv)-phiz*Zs-2*Ss/CC-Rs/CC-mu*eps_r,0,-eps*phi*Ss/Qv, mur*eps_r-Ss/CC,-phiz*Ss],
                  [eps*phi*Vs/Qv, -1/eta-ds-phiz*Zs, eps*phi*Ss/Qv+epsr*phir*Rs/Qv, epsr*phir*Vs/Qv, -phiz*Is],
                  [-phi*Vs/Qp, beta*Qv/(eta*Qp), -phi*Ss/Qp-phir*Rs/Qp-m-2*m2*Vs, -phir*Vs/Qp, 0],
                  [mu*eps_r-Rs/CC, 0, -epsr*phir*Rs/Qv, mur-ds-epsr*phir*Vs/Qv-phiz*Zs-2*Rs/CC-Ss/CC-mur*eps_r, -phiz*Rs],
                  [eps_z*phiz*Zs, eps_z*phiz*Zs, 0, eps_z*phiz*Zs, eps_z*phiz*(Is+Ss+Rs)-dz-2*dz2*Zs]])
    eigs=np.linalg.eig(Jac)
    return(eigs)


def Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Zs, Ss, Rs):
    J=[[ (-1/eta-d-phiz*Zs), phi*Ss/Qv+phir*Rs/Qv],[beta*Qv/(eta*Qp), -phi*Ss/Qp-phir*Rs/Qp-m]]
    eigs, e_vec=np.linalg.eig(J)
    eigs=[ei.real for ei in eigs]
    gr=np.where(np.absolute(eigs)==np.max(np.absolute(eigs)))
    gr=eigs[gr[0][0]]
    return gr

def MCT_coexistence_effects_IRS_on_VZ(S_aX, I_aX, R_aX,V_aX,Z_aX, S_c_aX, I_c_aX, R_c_aX, i1, i2,d, phiz, eps_z, dz,dz2, m, phi, Qp,Qv,beta, lat_per, phir, tinv):
    
    gZ0_aX=[]
    gZ0_cSI_aX=[]
    gZ0_cSR_aX=[]
    gZ0_cIR_aX=[]
    gZ0_cS_aX=[]
    gZ0_cI_aX=[]
    gZ0_cR_aX=[]
    gZ_aX=[]

    gV0_aX=[]
    gV0_cSI_aX=[]
    gV0_cSR_aX=[]
    gV0_cIR_aX=[]
    gV0_cS_aX=[]
    gV0_cI_aX=[]
    gV0_cR_aX=[]
    gV_aX=[]

    nsteps=len(S_aX)
    res=find_peaks(S_aX)
    res=res[0]
    if len(res)>1:
        p1=res[0]
        p2=res[1]
#        print('peaks found:')
#        print(p1)
#        print(p2)
#        print(S_aX[p1])
#        print(S_aX[p2])
    else:
        p1=0
        p2=len(S_aX)
    S_c_aX=np.mean(S_aX[p1:p2])
    I_c_aX=np.mean(I_aX[p1:p2])
    R_c_aX=np.mean(R_aX[p1:p2])
    for t in range(p1, p2):
    #for t in range(i2-i1):
        t_rand=random.randint(p1 , p2-1)
        
        if tinv=='zoop':
            gZ0_aX.append(phiz*eps_z*(S_c_aX+I_c_aX+R_c_aX)-dz)
            gZ0_cSI_aX.append(phiz*eps_z*(S_c_aX+I_c_aX+R_aX[t])-dz)
            gZ0_cSR_aX.append(phiz*eps_z*(S_c_aX+I_aX[t]+R_c_aX)-dz)
            gZ0_cIR_aX.append(phiz*eps_z*(S_aX[t]+I_c_aX+R_c_aX)-dz)
            gZ0_cS_aX.append(phiz*eps_z*(S_c_aX+I_aX[t]+R_aX[t])-dz)
            gZ0_cI_aX.append(phiz*eps_z*(S_aX[t]+I_c_aX+R_aX[t])-dz)
            gZ0_cR_aX.append(phiz*eps_z*(S_aX[t]+I_aX[t]+R_c_aX)-dz)
            gZ_aX.append(phiz*eps_z*(S_aX[t]+I_aX[t]+R_aX[t])-dz)
        else:
            gZ0_aX.append(phiz*eps_z*(S_c_aX+I_c_aX+R_c_aX)-dz-dz2*Z_aX[t])
            gZ0_cSI_aX.append(phiz*eps_z*(S_c_aX+I_c_aX+R_aX[t])-dz-dz2*Z_aX[t])
            gZ0_cSR_aX.append(phiz*eps_z*(S_c_aX+I_aX[t]+R_c_aX)-dz-dz2*Z_aX[t])
            gZ0_cIR_aX.append(phiz*eps_z*(S_aX[t]+I_c_aX+R_c_aX)-dz-dz2*Z_aX[t])
            gZ0_cS_aX.append(phiz*eps_z*(S_c_aX+I_aX[t]+R_aX[t])-dz-dz2*Z_aX[t])
            gZ0_cI_aX.append(phiz*eps_z*(S_aX[t]+I_c_aX+R_aX[t])-dz-dz2*Z_aX[t])
            gZ0_cR_aX.append(phiz*eps_z*(S_aX[t]+I_aX[t]+R_c_aX)-dz-dz2*Z_aX[t])
            gZ_aX.append(phiz*eps_z*(S_aX[t]+I_aX[t]+R_aX[t])-dz-dz2*Z_aX[t])

        #if tinv=='zoop':
#            print(beta*Qv/(lat_per*Qp)*I_c_aX/V_aX[t])
 #           print(-(m+phi*S_c_aX*1/Qp+phir*R_c_aX*1/Qp))
        gr0_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_c_aX, R_c_aX)
        gV0_aX.append(gr0_aX)
        gr0_cSI_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_c_aX, R_aX[t])
        gV0_cSI_aX.append(gr0_cSI_aX)
        gr0_cSR_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_c_aX, R_c_aX)
        gV0_cSR_aX.append(gr0_cSR_aX)
        gr0_cIR_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_aX[t], R_c_aX)
        gV0_cIR_aX.append(gr0_cIR_aX)
        gr0_cS_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_c_aX, R_aX[t])
        gV0_cS_aX.append(gr0_cS_aX)
        gr0_cI_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_aX[t], R_aX[t])
        gV0_cI_aX.append(gr0_cI_aX)
        gr0_cR_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_aX[t], R_c_aX)
        gV0_cR_aX.append(gr0_cR_aX)
        gr_aX=Jacob_IV(phi,phir, Qv,Qp, lat_per, d,m, phiz, beta, Z_aX[t], S_aX[t], R_aX[t])
        gV_aX.append(gr_aX)
        #else:
        #    gV0_aX.append(-(m+phi*S_c_aX*1/Qp+phir*R_c_aX*1/Qp))
        #    gV0_cSI_aX.append(-(m+phi*S_c_aX*1/Qp+phir*R_aX[t]*1/Qp))
        #    gV0_cSR_aX.append(-(m+phi*S_c_aX*1/Qp+phir*R_c_aX*1/Qp))
        #    gV0_cIR_aX.append(-(m+phi*S_aX[t]*1/Qp+phir*R_c_aX*1/Qp))
        #    gV0_cS_aX.append(-(m+phi*S_c_aX*1/Qp+phir*R_aX[t]*1/Qp))
        #    gV0_cI_aX.append(-(m+phi*S_aX[t]*1/Qp+phir*R_aX[t]*1/Qp))
        #    gV0_cR_aX.append(-(m+phi*S_aX[t]*1/Qp+phir*R_c_aX*1/Qp))
        #    gV_aX.append(-(m+phi*S_aX[t]*1/Qp+phir*R_aX[t]*1/Qp))

    epsZ_0_aX=np.mean(gZ0_aX)
    epsZ_S_aX=np.mean(gZ0_cIR_aX)-epsZ_0_aX
    epsZ_I_aX=np.mean(gZ0_cSR_aX)-epsZ_0_aX
    epsZ_R_aX=np.mean(gZ0_cSI_aX)-epsZ_0_aX
    epsZ_SI_aX= np.mean(gZ0_cR_aX)-(epsZ_0_aX+epsZ_S_aX+epsZ_I_aX)
    epsZ_SR_aX= np.mean(gZ0_cI_aX)-(epsZ_0_aX+epsZ_S_aX+epsZ_R_aX)
    epsZ_IR_aX= np.mean(gZ0_cS_aX)-(epsZ_0_aX+epsZ_I_aX+epsZ_R_aX)
    epsZ_SIR_aX= np.mean(gZ_aX)-(epsZ_0_aX+epsZ_S_aX+epsZ_I_aX+epsZ_R_aX+epsZ_SI_aX+epsZ_SR_aX+epsZ_IR_aX)

    epsV_0_aX=np.mean(gV0_aX)
    epsV_S_aX=np.mean(gV0_cIR_aX)-epsV_0_aX
    epsV_I_aX=np.mean(gV0_cSR_aX)-epsV_0_aX
    epsV_R_aX=np.mean(gV0_cSI_aX)-epsV_0_aX
    epsV_SI_aX= np.mean(gV0_cR_aX)-(epsV_0_aX+epsV_S_aX+epsV_I_aX)
    epsV_SR_aX= np.mean(gV0_cI_aX)-(epsV_0_aX+epsV_S_aX+epsV_R_aX)
    epsV_IR_aX= np.mean(gV0_cS_aX)-(epsV_0_aX+epsV_I_aX+epsV_R_aX)
    epsV_SIR_aX= np.mean(gV_aX)-(epsV_0_aX+epsV_S_aX+epsV_I_aX+epsV_R_aX+epsV_SI_aX+epsV_SR_aX+epsV_IR_aX)
    
    effs=[epsZ_0_aX, epsZ_S_aX, epsZ_I_aX,epsZ_R_aX, epsZ_SI_aX,epsZ_SR_aX,epsZ_IR_aX,epsZ_SIR_aX,
         epsV_0_aX, epsV_S_aX, epsV_I_aX,epsV_R_aX, epsV_SI_aX,epsV_SR_aX,epsV_IR_aX,epsV_SIR_aX]
    return effs
    
    
def MCT_analysis_SIVRZ(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,eps,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,dt, ndays, ntot):
    # absent/invading zoop
    init_conditions=[0.001, 0.001,0.001,0.001,0]
    alpha=1
    result_az=simulation_SIVRZ_rk4(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha, dt, ndays, init_conditions)
    
    Nt=len(result_az[0])
    i1=Nt-round(365/dt)
    i2=Nt
    S_az=result_az[0][i1:i2]
    I_az=result_az[1][i1:i2]
    V_az=result_az[2][i1:i2]
    R_az=result_az[3][i1:i2]
    Z_az=result_az[4][i1:i2]

    S_c_az=np.mean(S_az)
    I_c_az=np.mean(I_az)
    R_c_az=np.mean(R_az)
    
    tinv='zoop'
    effects_az=MCT_coexistence_effects_IRS_on_VZ(S_az, I_az,R_az, V_az,Z_az,S_c_az, I_c_az, R_c_az, i1, i2,ds, phiz, eps_z, dz,dz2, m, phi, Qp,Qv,beta, eta, phir, tinv)
    
    # absent/invading virus
    init_conditions=[0.001, 0,0,0.001,0.001]
    result_av=simulation_SIVRZ_rk4(mu, mui, eta, beta, phi, ds, m,m2, Qv, Qp,  eps,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,alpha, dt, ndays, init_conditions)

    S_av=result_av[0][i1:i2]
    I_av=result_av[1][i1:i2]
    V_av=result_av[2][i1:i2]
    R_av=result_av[3][i1:i2]
    Z_av=result_av[4][i1:i2]

    S_c_av=np.mean(S_av)
    I_c_av=np.mean(I_av)
    R_c_av=np.mean(R_av)
    
    tinv='virus'
    effects_av=MCT_coexistence_effects_IRS_on_VZ(S_av, I_av,R_av,V_av,Z_av,S_c_av, I_c_av, R_c_av,i1, i2,ds, phiz, eps_z, dz,dz2, m, phi, Qp,Qv,beta, eta, phir, tinv)
    
    delta_effects=[]
    for i in range(ntot-1):
        de=effects_av[i+ntot-1]-effects_av[i]
        delta_effects.append(de)
    for i in range(ntot-1):
        de0=effects_az[i]-effects_az[i+ntot-1]
        delta_effects.append(de0)
    return delta_effects

def coexistence_analysis_SIVRZ(S, I, V,R, Z,mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,eps,epso,eps_r,eps_lr,phir,mur,phiz,eps_z, dz,dz2,CC,dt, ndays, ntot):
    gZ0=[]
    gZ0_cSI=[]
    gZ0_cSR=[]
    gZ0_cIR=[]
    gZ0_cS=[]
    gZ0_cI=[]
    gZ0_cR=[]
    gZ=[]

    gV0=[]
    gV0_cSI=[]
    gV0_cSR=[]
    gV0_cIR=[]
    gV0_cS=[]
    gV0_cI=[]
    gV0_cR=[]
    gV=[]

    lp=eta
    nsteps=len(S)
    res=find_peaks(S)
    res=res[0]
    if len(res)>1:
        p1=res[0]
        p2=res[1]
        #print('peaks found:')
        #print(p1)
        #print(p2)
        #print(S[p1])
        #print(S[p2])
    else:
        p1=0
        p2=len(S)
    S_c=np.mean(S[p1:p2])
    I_c=np.mean(I[p1:p2])
    R_c=np.mean(R[p1:p2])
    for t in range(p1, p2):
    #for t in range(i2-i1):
        t_rand=random.randint(p1 , p2-1)

        gZ0.append(phiz*eps_z*(S_c+I_c+R_c)-dz-dz2*Z[t])
        gZ0_cSI.append(phiz*eps_z*(S_c+I_c+R[t])-dz-dz2*Z[t])
        gZ0_cSR.append(phiz*eps_z*(S_c+I[t]+R_c)-dz-dz2*Z[t])
        gZ0_cIR.append(phiz*eps_z*(S[t]+I_c+R_c)-dz-dz2*Z[t])
        gZ0_cS.append(phiz*eps_z*(S_c+I[t]+R[t])-dz-dz2*Z[t])
        gZ0_cI.append(phiz*eps_z*(S[t]+I_c+R[t])-dz-dz2*Z[t])
        gZ0_cR.append(phiz*eps_z*(S[t]+I[t]+R_c)-dz-dz2*Z[t])
        gZ.append(phiz*eps_z*(S[t]+I[t]+R[t])-dz-dz2*Z[t])

        #if tinv=='zoop':
#            print(beta*Qv/(lat_per*Qp)*I_c_aX/V_aX[t])
 #           print(-(m+phi*S_c_aX*1/Qp+phir*R_c_aX*1/Qp))
        gr0=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S_c/Qp-phir*R_c/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S_c, R_c)
        gV0.append(gr0)
        gr0_cSI=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S_c/Qp-phir*R[t]/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S_c, R[t])
        gV0_cSI.append(gr0_cSI)
        gr0_cSR=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S_c/Qp-phir*R_c/Qp-m2*V[t]   #Jacob_IV(phi,phir, Qv,Qp,eta, d,m, phiz, beta, Z[t], S_c, R_c)
        gV0_cSR.append(gr0_cSR)
        gr0_cIR=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S[t]/Qp-phir*R_c/Qp-m2*V[t] #Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S[t], R_c)
        gV0_cIR.append(gr0_cIR)
        gr0_cS=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S_c/Qp-phir*R[t]/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S_c, R[t])
        gV0_cS.append(gr0_cS)
        gr0_cI=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S[t]/Qp-phir*R[t]/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S[t], R[t])
        gV0_cI.append(gr0_cI)
        gr0_cR=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S[t]/Qp-phir*R_c/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S[t], R_c)
        gV0_cR.append(gr0_cR)
        gr=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S[t]/Qp-phir*R[t]/Qp-m2*V[t]#Jacob_IV(phi,phir, Qv,Qp, eta, d,m, phiz, beta, Z[t], S[t], R[t])
        gV.append(gr)

    epsZ_0=np.mean(gZ0)
    epsZ_S=np.mean(gZ0_cIR)-epsZ_0
    epsZ_I=np.mean(gZ0_cSR)-epsZ_0
    epsZ_R=np.mean(gZ0_cSI)-epsZ_0
    epsZ_SI= np.mean(gZ0_cR)-(epsZ_0+epsZ_S+epsZ_I)
    epsZ_SR= np.mean(gZ0_cI)-(epsZ_0+epsZ_S+epsZ_R)
    epsZ_IR= np.mean(gZ0_cS)-(epsZ_0+epsZ_I+epsZ_R)
    epsZ_SIR= np.mean(gZ)-(epsZ_0+epsZ_S+epsZ_I+epsZ_R+epsZ_SI+epsZ_SR+epsZ_IR)

    epsV_0=np.mean(gV0)
    epsV_S=np.mean(gV0_cIR)-epsV_0
    epsV_I=np.mean(gV0_cSR)-epsV_0
    epsV_R=np.mean(gV0_cSI)-epsV_0
    epsV_SI= np.mean(gV0_cR)-(epsV_0+epsV_S+epsV_I)
    epsV_SR= np.mean(gV0_cI)-(epsV_0+epsV_S+epsV_R)
    epsV_IR= np.mean(gV0_cS)-(epsV_0+epsV_I+epsV_R)
    epsV_SIR= np.mean(gV)-(epsV_0+epsV_S+epsV_I+epsV_R+epsV_SI+epsV_SR+epsV_IR)

    effs=[epsZ_0, epsZ_S, epsZ_I,epsZ_R, epsZ_SI,epsZ_SR,epsZ_IR,epsZ_SIR,
         epsV_0, epsV_S, epsV_I,epsV_R, epsV_SI,epsV_SR,epsV_IR,epsV_SIR]
    return effs
