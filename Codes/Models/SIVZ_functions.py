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
from generic_functions import *

# integartion
def dSIVZ(St, It,Vt, Zt, dt,mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph):

    dS=dt*( (mu-d-eps*phi*Vt/Qv-St/CC)*St +epso*It-phiz*Zt*St)
    dI=dt*((mui-d-(1/eta))*It+eps*phi*Vt*St/Qv-epso*It-alph*phiz*Zt*It)
    dV=dt*((1/eta)*beta*It*(Qv/Qp)-(m+phi*St*1/Qp)*Vt-m2*Vt*Vt)
    dZ=dt*((phiz*eps_z*(St+It)-dz-dz2*Zt)*Zt)
    mZ=phiz*Zt
    if It+St>0:
        pi=It/(It+St)
        mV=pi/eta
    else:
        mV=0
    if Zt !=0:
        gZ=phiz*eps_z*(St+It)-dz
    else:
        gZ=0
    if Vt !=0:
        gV=(1/eta)*beta*It/Vt*(Qv/Qp)-(m+phi*St*1/Qp)-m2*Vt
    else:
        gV=0
    npp=(mu-d)*St
    return dS, dI, dV, dZ, mZ, mV, gZ, gV, npp

# runge lutta 4
def rk4_SIVZ(St, It,Vt, Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph):
    k1S, k1I, k1V, k1Z, m1Z, m1V, g1Z, g1V, npp1=dSIVZ(St, It,Vt, Zt,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph)
    k2S, k2I, k2V, k2Z, m2Z, m2V, g2Z, g2V, npp2=dSIVZ(St+0.5*k1S, It+0.5*k1I,Vt+0.5*k1V, Zt+0.5*k1Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph)
    k3S, k3I, k3V, k3Z, m3Z, m3V, g3Z, g3V, npp3=dSIVZ(St+0.5*k2S, It+0.5*k2I,Vt+0.5*k2V, Zt+0.5*k2Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph)
    k4S, k4I, k4V, k4Z, m4Z, m4V, g4Z, g4V, npp4=dSIVZ(St+k3S, It+k3I,Vt+k3V, Zt+k3Z,dt, mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph)
    Stn=St+(1.0/6.0)*(k1S + 2*k2S + 2*k3S + k4S)
    Itn=It+(1.0/6.0)*(k1I + 2*k2I + 2*k3I + k4I)
    Vtn=Vt+(1.0/6.0)*(k1V + 2*k2V + 2*k3V + k4V)
    Ztn=Zt+(1.0/6.0)*(k1Z + 2*k2Z + 2*k3Z + k4Z)
    mZ=(1.0/6.0)*(m1Z+2*m2Z+2*m3Z+m4Z)
    mV=(1.0/6.0)*(m1V+2*m2V+2*m3V+m4V)
    gZ=(1.0/6.0)*(g1Z+2*g2Z+2*g3Z+g4Z)
    gV=(1.0/6.0)*(g1V+2*g2V+2*g3V+g4V)
    npp=(1.0/6.0)*(npp1+2*npp2+2*npp3+npp4)
    return Stn, Itn, Vtn, Ztn, mZ, mV, gZ, gV, npp

# simulation
def simulation_SIVZ_rk4(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  
                        eps,epso,phiz,eps_z, dz,dz2,CC,alph, dt, ndays, init_conditions):
    #inittial conditions
    S0=init_conditions[0]
    I0=init_conditions[1]
    V0=init_conditions[2]
    Z0=init_conditions[3]
    S=[S0]
    I=[I0]
    V=[V0]
    Z=[Z0]
    mVs=[0]
    mZs=[0]
    gVs=[0]
    gZs=[0]
    npps=[0]
    time_sim=range(round(ndays/dt))


    for t in time_sim:
        #print(t)
        St=S[t]
        It=I[t]
        Vt=V[t]
        Zt=Z[t]
        Stn, Itn, Vtn, Ztn, mZ, mV, gZ, gV, npp=rk4_SIVZ(St, It, Vt, Zt, dt, mu, mui, eta,
                                                    beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, alph)
        Stn=max([0,Stn])
        Itn=max([0,Itn])
        Vtn=max([0,Vtn])
        Ztn=max([0,Ztn])

        mZs.append(mZ)
        mVs.append(mV)
        gZs.append(gZ)
        gVs.append(gV)
        npps.append(npp)

        S.append(Stn)
        I.append(Itn)
        V.append(Vtn)
        Z.append(Ztn)
    return S, I, V, Z, mZs, mVs, gZs, gVs, npps

# equilibrium in the case dz2!=0 and dv2=0
def equilibrium_SIVZ(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC):
    g=phiz
    eta=lp
    a=mu-d+g*(dz/dz2)
    b=1/CC+g*g*eps_z/dz2
    T=(phi*beta)/(m*eta*Qp)+g*g*eps_z/dz2
    S_star=(1/eta+d+g*g*eps_z/dz2*(a/T-dz))/((phi*beta)/(m*eta*Qp)+g*g*eps_z/dz2*(b/T-1))
    I_star=(a-b*S_star)/T

    V_star=beta*Qv*I_star/((m+phi*S_star/Qp)*eta*Qp)
    Z_star=(eps_z*g*(I_star+S_star)-dz)/dz2
    
    return S_star, I_star, V_star, Z_star

# equilibrium in the case dz2!=0 and dv2!=0
def equilibrium_SIVZ_m2(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,  eps, epso,phiz,eps_z,dz,dz2, CC, n_update):
    # first approximation
    g=phiz
    a0=mu-d+g*(dz/dz2)+phi*eps/Qp*(m/m2)
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
        for i in range(n_update):
            S_star0=S_star_SVZ
            Z_star0=(eps_z*g*(I_star0+S_star0)-dz)/dz2
            V_star0=mt.sqrt((m+phi*S_star0/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star0)/(2*m2)-m/(2*m2)-phi*S_star0/(Qp*2*m2)

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
                # exclude the point if I_star is found to be too low 
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

    # decide which approximation to use
    if ((v2>10*v1) and I_star0>0): #or phi/Qp<0.5:
        case=1
    else:
        case=2
    if case==1 and I_star0>0:
        A=beta*Qv/(lp*Qp*m2)

        denom=1/CC+g*g*eps_z/dz2

        B=(phi/Qv*mt.sqrt(A))/denom
        C=(mu-d+g*dz/dz2)/denom
        D=((g**2)*eps_z/dz2)/denom
        E=-d-1/lp+g*dz/dz2
        F=eps*phi*phi/(Qv*Qp*2*m2)
        
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
            V_star=mt.sqrt((m+phi*S_star/Qp)**2+4*m2*beta*Qv/(lp*Qp)*I_star)/(2*m2)-m/(2*m2)-phi*S_star/(Qp*2*m2) 
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
    elif case==2 and lp<=7:
        S_star=S_star0
        I_star=I_star0
        V_star=V_star0
        Z_star=Z_star0
    elif case==2 and lp>7:
        S_star=np.nan
        I_star=np.nan
        V_star=np.nan
        Z_star=np.nan

    return S_star, I_star, V_star, Z_star, case

# theoretical equilibrium in the case dz2!=0 and dv2=0 and alpha=0
def equilibrium_SIVZ_alpha(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, alpha):
    g=phiz
    eta=lp
    a=mu-d+g*(dz/dz2)
    b=1/CC+g*g*eps_z/dz2
    T=(phi*beta)/(m*eta*Qp)+g*g*eps_z*alpha/dz2
    try:
        S_star=(1/eta+d+g*g*eps_z/dz2*(a*alpha/T-dz))/((phi*beta)/(m*eta*Qp)+g*g*eps_z/dz2*(b*alpha/T-1))
    except ZeroDivisionError:
        S_star=0
    I_star=(a-b*S_star)/T
    V_star=beta*Qv*I_star/(m*eta*Qp)
    Z_star=(eps_z*g*(alpha*I_star+S_star)-dz)/dz2
    return S_star, I_star, V_star, Z_star

# theoretical equilibrium in the case dz2!=0 and dv2=0 and loss of infection !=0 (not included in the paper)
def equilibrium_SIVZ_LOI(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, eff):
    g=phiz
    gamma=float(eff)
    eta=lp

    a=1/eta+d+gamma-g*dz/dz2
    b=g*g*eps_z/dz2
    T=(beta*phi)/(m*eta*Qp)-g*g*eps_z/dz2

    if T>0:
        c=mu-d+g*dz/dz2
        d0=g*g*eps_z/dz2+1/CC
        e=g*g*eps_z/dz2+(beta*phi)/(m*eta*Qp)

        A = -d0*pow(b/T,2)-b*e/T
        B = (b*c-a*e)/T-2*a*b*d0/pow(T,2)+gamma
        C = -d0*pow(a/T,2)+a*c/T
        if pow(B,2)-4*A*C>0:
            sol1 = (-B+mt.sqrt(pow(B,2)-4*A*C))/(2*A)
            sol2=  (-B-mt.sqrt(pow(B,2)-4*A*C))/(2*A)
            I_star = max([sol1, sol2])
            S_star = (a+b*I_star)/T
            V_star_ap=beta*Qv*I_star/(m*eta*Qp)
            V_star=beta*Qv*I_star/(eta*Qp)/(m+phi*S_star/Qp)
            Z_star=(eps_z*g*(I_star+S_star)-dz)/dz2
        else:
            I_star =0
            S_star =0
            V_star=0
            Z_star=0
    else:
        I_star =0
        S_star =0
        V_star=0
        Z_star=0
    return S_star, I_star, V_star, Z_star

# theoretical equilibrium in the case dz2!=0 and dv2=0 and loss of infection !=0 and alpha!=0 (not included in the paper)
def equilibrium_SIVZ_LOI_alpha(mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC, eff, alpha):
   gamma=float(eff)
   eta=lp
   g=phiz

   a=1/eta+d+gamma-g*dz/dz2
   b=g*g*eps_z*alpha/dz2
   T=(beta*phi)/(m*eta*Qp)-g*g*eps_z/dz2
   if T>0:
       c=mu-d+g*dz/dz2
       d0=g*g*eps_z/dz2+1/CC
       e=g*g*eps_z*alpha/dz2+(beta*phi)/(m*eta*Qp)

       A = -d0*pow(b/T,2)-b*e/T
       B = (b*c-a*e)/T-2*a*b*d0/pow(T,2)+gamma
       C = -d0*pow(a/T,2)+a*c/T
       if pow(B,2)-4*A*C>0 and A!=0:
           sol1 = (-B+mt.sqrt(pow(B,2)-4*A*C))/(2*A)
           sol2=  (-B-mt.sqrt(pow(B,2)-4*A*C))/(2*A)
           I_star = max([sol1, sol2])
           S_star = (a+b*I_star)/T
           V_star_ap=beta*Qv*I_star/(m*eta*Qp)
           V_star=beta*Qv*I_star/(eta*Qp)/(m+phi*S_star/Qp)
           Z_star=(eps_z*g*(alpha*I_star+S_star)-dz)/dz2
       else:
           I_star = -C/B
           S_star = a/T
           V_star_ap=beta*Qv*I_star/(m*eta*Qp)
           V_star=beta*Qv*I_star/(eta*Qp)/(m+phi*S_star/Qp)
           Z_star=(eps_z*g*(alpha*I_star+S_star)-dz)/dz2
   else:
       I_star=0
       S_star=0
       V_star=0
       Z_star=0
   return S_star, I_star, V_star, Z_star

# check stable state in the case dz2!=0 and dv2!=0 (m2 in the code)
def check_stable_states_SIVZ_m2(mu, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, n_state):
    epsilon=0
    lp=eta
    g=phiz
    possible_init_values=[]
    #S, I, V, Z
    S_star, I_star, V_star, Z_star,case=equilibrium_SIVZ_m2(mu, 0, lp, beta, phi, d, m,m2, Qv, Qp,  eps, epso,phiz,eps_z,dz,dz2, CC, 1000)
    ast=np.nan
    ost=np.nan
    fst=np.nan
    if V_star>0 and Z_star>0 and S_star>0 and I_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
        ast, ost, fst=check_eigs(eigs, eig_vecs, ast, ost)
    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # S, I, V, 0
    S_star, I_star, V_star, Z_star,case=equilibrium_SIVZ_m2(mu, 0, lp, beta, phi, d, m,m2, Qv, Qp,  eps, epso,0,0,0,dz2, CC, 1000) # we keep dz2 as it appears in denominators of terms related to Z
    Z_star=epsilon
    ast1=np.nan
    ost1=np.nan
    fst1=np.nan
    if  S_star>0 and V_star>0 and I_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star,Z_star)
        ast1, ost1, fst1=check_eigs(eigs, eig_vecs, ast1, ost1)  
    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # S, 0 , 0, Z
    S_star=(mu-d+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*S_star-dz)/dz2
    I_star=epsilon
    V_star=epsilon
    ast2=np.nan
    ost2=np.nan
    fst2=np.nan
    if  Z_star>0 and S_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
        ast2, ost2, fst2=check_eigs(eigs, eig_vecs, ast2, ost2)


    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # S, 0,0,0
    S_star=(mu-d)*CC
    Z_star=epsilon
    I_star=epsilon
    V_star=epsilon
    ast3=np.nan
    ost3=np.nan
    fst3=np.nan
    if  S_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star,Z_star)
        ast3, ost3, fst3=check_eigs(eigs, eig_vecs, ast3, ost3)


    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # 0, 0, 0, 0
    ast4=1
    st4=1
    ost4=1
    fst4=0
    S_star=epsilon
    I_star=epsilon
    V_star=epsilon
    Z_star=epsilon
    eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC,S_star, I_star, V_star, Z_star)
    possible_init_values.append([epsilon, epsilon, epsilon, epsilon])
    ast4, ost4, fst4=check_eigs(eigs, eig_vecs, ast4, ost4)

    a_stable=[ast, ast1,ast2,ast3,ast4]
    o_stable=[ost, ost1,ost2,ost3,ost4]
    f_stable=[fst, fst1,fst2,fst3,fst4]

    init_conditions=[0.001, 0.001,0.001,0.001]
    for l in range(n_state):
        if a_stable[l]==1:
            init_conditions=possible_init_values[l]
            break

    return a_stable, o_stable,f_stable,init_conditions

# check stable state in the case dz2!=0 and dv2=0 (m2 in the code)
def check_stable_states_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps,epso,phiz,eps_z, dz,dz2,CC, n_state):
    possible_init_values=[]
    mui=0
    epsilon=0
    if epso==0:
        # full coexistance S, I, V, Z
        g=phiz
        ast=np.nan
        ost=np.nan
        fst=np.nan
        S_star, I_star, V_star, Z_star=equilibrium_SIVZ(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC)
        if V_star>0 and Z_star>0 and S_star>0 and I_star>0:
            eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
            ast, ost, fst=check_eigs(eigs, eig_vecs, ast, ost)

        possible_init_values.append([S_star, I_star, V_star, Z_star])

        # S, I, V, 0
        ast1=np.nan
        ost1=np.nan
        fst1=np.nan
        S_star0, I_star0, V_star0, Z_star0=equilibrium_SIVZ(mu, mui, eta, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,0,0,0,dz2, CC)
        if  S_star0>0 and V_star0>0 and I_star0>0:
            eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star0, I_star0, V_star0,Z_star0)
            ast1, ost1, fst1=check_eigs(eigs, eig_vecs, ast1, ost1)


        possible_init_values.append([S_star0, I_star0, V_star0, Z_star0])
    elif epso!=0:
        # full coexistance S, I, V, Z
        g=phiz
        gamma=float(epso)

        a=1/eta+d+gamma-g*dz/dz2
        b=g*g*eps_z/dz2
        T=(beta*phi)/(m*eta*Qp)-g*g*eps_z/dz2
        ast=np.nan
        ost=np.nan
        fst=np.nan
        S_star=0
        I_star=0
        V_star=0
        Z_star=0
        if T>0:
            c=mu-d+g*dz/dz2
            d0=g*g*eps_z/dz2+1/CC
            e=g*g*eps_z/dz2+(beta*phi)/(m*eta*Qp)

            A = -d0*pow(b/T,2)-b*e/T
            B = (b*c-a*e)/T-2*a*b*d0/pow(T,2)+gamma
            C = -d0*pow(a/T,2)+a*c/T
            if pow(B,2)-4*A*C>0:
                sol1 = (-B+mt.sqrt(pow(B,2)-4*A*C))/(2*A)
                sol2=  (-B-mt.sqrt(pow(B,2)-4*A*C))/(2*A)
                I_star = max([sol1, sol2])
                S_star = (a+b*I_star)/T
                V_star_ap=beta*Qv*I_star/(m*eta*Qp)
                V_star=beta*Qv*I_star/(eta*Qp)/(m+phi*S_star/Qp)
                Z_star=(eps_z*g*(I_star+S_star)-dz)/dz2
                ast=np.nan
                st=np.nan
                ost=np.nan
                if V_star>0 and Z_star>0 and S_star>0:
                    eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
                    ast, ost, fst=check_eigs(eigs, eig_vecs, ast, ost)


        possible_init_values.append([S_star, I_star, V_star, Z_star])

        # S, I, V, 0
        S_star0=(d+1/eta+gamma)/(phi*beta/(m*eta*Qp))
        I_star0=(mu-d-S_star0/CC)/(phi*beta/(m*eta*Qp)-gamma/S_star0)
        V_star_ap=beta*Qv*I_star0/(m*eta*Qp)
        V_star0=beta*Qv*I_star0/(eta*Qp)/(m+phi*S_star0/Qp)
        Z_star0=epsilon
        ast1=np.nan
        st1=np.nan
        ost1=np.nan
        fst1=np.nan
        if  S_star0>0 and V_star0>0 and I_star0>0:
            eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star0, I_star0, V_star0,Z_star0)
            ast1, ost1, fst1=check_eigs(eigs, eig_vecs, ast1, ost1)


        possible_init_values.append([S_star0, I_star0, V_star0, Z_star0])

    # S, 0 , 0, Z
    S_star=(mu-d+g*dz/dz2)/(1/CC+g*g*eps_z/dz2)
    Z_star=(eps_z*g*S_star-dz)/dz2
    I_star=epsilon
    V_star=epsilon
    ast2=np.nan
    ost2=np.nan
    fst2=np.nan
    if  Z_star>0 and S_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star, Z_star)
        ast2, ost2, fst2=check_eigs(eigs, eig_vecs, ast2, ost2)


    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # S, 0,0,0
    S_star=(mu-d)*CC
    Z_star=epsilon
    I_star=epsilon
    V_star=epsilon
    ast3=np.nan
    ost3=np.nan
    fst3=np.nan
    if  S_star>0:
        eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC, S_star, I_star, V_star,Z_star)
        ast3, ost3, fst3=check_eigs(eigs, eig_vecs, ast3, ost3)


    possible_init_values.append([S_star, I_star, V_star, Z_star])

    # 0, 0, 0, 0
    ast4=1
    st4=1
    ost4=1
    fst4=0
    S_star=epsilon
    I_star=epsilon
    V_star=epsilon
    Z_star=epsilon
    eigs, eig_vecs=Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,epso,phiz,eps_z, dz,dz2,CC,S_star, I_star, V_star, Z_star)
    possible_init_values.append([epsilon, epsilon, epsilon, epsilon])
    ast4, ost4, fst4=check_eigs(eigs, eig_vecs, ast4, ost4)

    a_stable=[ast, ast1,ast2,ast3,ast4]
    o_stable=[ost, ost1,ost2,ost3,ost4]
    f_stable=[fst, fst1,fst2,fst3,fst4]

    init_conditions=[0.001, 0.001,0.001,0.001]
    for l in range(n_state):
        if a_stable[l]==1:
            init_conditions=possible_init_values[l]
            break

    return a_stable, o_stable,f_stable,init_conditions

# time series plot 1
def make_1_plot_SIVZ(S, I, V, Z,i1,i2, title,dt, pp, Qp, Qv, Qz, plot_peaks):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    time_sim=range(i1, i2)
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    bottom, top = plt.ylim()
    plt.yscale('log')
    leg_vec=['Susceptible', 'Infected', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    if plot_peaks=='yes':
        V_peak, V_peak_l, V_peak_b=find_low_and_high(V[i1:i2])
        Z_peak, Z_peak_l, Z_peak_b=find_low_and_high(Z[i1:i2])
        plt.axvline(x=V_peak*dt+time_sim[0]*dt, color='#1f77b4', linestyle='-', linewidth=1)
        plt.axvline(x=V_peak_l*dt+time_sim[0]*dt, color='#1f77b4', linestyle='-', linewidth=1)
        plt.axvline(x=Z_peak*dt+time_sim[0]*dt, color='red', linestyle='-', linewidth=1)
        plt.axvline(x=Z_peak_l*dt+time_sim[0]*dt, color='red', linestyle='-', linewidth=1)
    pp.savefig()

# time series plot 2
def make_plots_SIVZ(S, I, V, Z,i1,i2, title,dt, pp, Qp, Qv, Qz):
    plt.rcParams['lines.linewidth'] = 3
    #time step in day
    time_sim=range(i1, i2)
    # molar time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, S[i1:i2], color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, I[i1:i2], color='#ff7f0e')
    #plt.plot(np.array(time_sim)*dt, R[0:len(R)-1])
    plt.plot(np.array(time_sim)*dt, V[i1:i2], color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, Z[i1:i2], color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (umolN.L-1)')
    bottom, top = plt.ylim()
    plt.ylim((1e-7,0.5))
    plt.yscale('log')
    leg_vec=['Susceptible', 'Infected', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()

    #count time series
    fig, ax = plt.subplots(figsize=(4,4))
    plt.plot(np.array(time_sim)*dt, np.array(S[i1:i2])/Qp, color='#2ca02c')
    plt.plot(np.array(time_sim)*dt, np.array(I[i1:i2])/Qp, color='#ff7f0e')
    plt.plot(np.array(time_sim)*dt, np.array(V[i1:i2])/Qv, color='#1f77b4')
    plt.plot(np.array(time_sim)*dt, np.array(Z[i1:i2])/Qz, color='red')
    plt.xlabel('Days')
    plt.ylabel('Concentration (ind.L-1)')
    plt.yscale('log')
    mini=1
    maxi=5e10
    plt.ylim((mini,maxi))
    leg_vec=['Susceptible', 'Infected', 'Virus', 'Zooplankton']
    plt.legend(leg_vec,bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0., ncol=4, fontsize=6)
    plt.title(title)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    pp.savefig()


# jacobian of the SIVZ model
def Jacobian_SIVZ(mu, eta, beta, phi, d, m,m2, Qv, Qp,  eps,
                       epso,phiz,eps_z, dz,dz2,CC, Ss, Is, Vs, Zs):
    Jac=np.array([[mu-d-phi*Vs/Qv-phiz*Zs-2*Ss/CC, epso, -phi*Ss/Qv,-phiz*Ss] ,
              [phi*Vs/Qv, -1/eta-d-phiz*Zs-epso, phi*Ss/Qv, -phiz*Is],
              [-phi*Vs/Qp, beta*Qv/(eta*Qp), -phi*Ss/Qp-m-2*m2*Vs, 0],
              [eps_z*phiz*Zs, eps_z*phiz*Zs, 0, eps_z*phiz*(Is+Ss)-dz-2*dz2*Zs]])
    eigs=np.linalg.eig(Jac)
    return(eigs)



# jacobian of the IV subsystem (growth rate of the virus)             
def Jacob_IV(Qv,Qp, eta, d,m, m2,phiz, epso, beta, phi,Vs, Zs, Is, Ss):
    J=[[ (-1/eta-d-phiz*Zs-epso), phi*Ss/Qv],[beta*Qv/(eta*Qp), -phi*Ss/Qp-m-2*m2*Vs]]
    eigs, e_vec=np.linalg.eig(J)
    eigs=[ei.real for ei in eigs]
    gr=max(eigs)
    return gr

# calculate effects of Modern Coexistence Theory form Ellner et al 2019
def MCT_coexistence_effects_IS_on_VZ(S_aX, I_aX,Z_aX,V_aX, S_c_aX, I_c_aX, i1, i2,mu,eps, epso, ds,phiz, eps_z, dz,dz2, m,m2, phi, Qp,Qv, beta, lat_per,alph,CC, tinv):
    gZ_aX=[]
    gZ0_aX=[]
    gZ0_cS_aX=[]
    gZ0_cI_aX=[]
    gZ0_ncovIS_aX=[]

    gV_aX=[]
    gV0_aX=[]
    gV0_cS_aX=[]
    gV0_cI_aX=[]
    gV0_ncovIS_aX=[]
    nsteps=len(S_aX)
    for t in range(nsteps):
        t_rand=random.randint(0 , nsteps-1)
        if tinv=='zoop':
            gZ_aX.append(phiz*eps_z*(S_aX[t]+alph*I_aX[t])-dz)
            gZ0_aX.append(phiz*eps_z*(S_c_aX+alph*I_c_aX)-dz)
            gZ0_cS_aX.append(phiz*eps_z*(S_c_aX+alph*I_aX[t])-dz)
            gZ0_cI_aX.append(phiz*eps_z*(S_aX[t]+alph*I_c_aX)-dz)
            gZ0_ncovIS_aX.append(phiz*eps_z*(S_aX[t]+alph*I_aX[t_rand])-dz)
        else:
            gZ_aX.append(phiz*eps_z*(S_aX[t]+alph*I_aX[t])-dz-dz2*Z_aX[t])
            gZ0_aX.append(phiz*eps_z*(S_c_aX+alph*I_c_aX)-dz-dz2*Z_aX[t])
            gZ0_cS_aX.append(phiz*eps_z*(S_c_aX+alph*I_aX[t])-dz-dz2*Z_aX[t])
            gZ0_cI_aX.append(phiz*eps_z*(S_aX[t]+alph*I_c_aX)-dz-dz2*Z_aX[t])
            gZ0_ncovIS_aX.append(phiz*eps_z*(S_aX[t]+alph*I_aX[t_rand])-dz-dz2*Z_aX[t])
        
        gr_aX=Jacob_IV(Qv,Qp, lat_per, ds,m,m2, phiz, epso,beta,phi, V_aX[t], Z_aX[t], I_aX[t], S_aX[t])
        gV_aX.append(gr_aX)
        
        gr0_aX=Jacob_IV(Qv,Qp, lat_per, ds,m,m2,  phiz, epso,beta,phi, V_aX[t], Z_aX[t], I_c_aX, S_c_aX)
        gV0_aX.append(gr0_aX)

        gr0_cS_aX=Jacob_IV(Qv,Qp, lat_per, ds,m,m2,  phiz, epso,beta,phi, V_aX[t], Z_aX[t], I_aX[t], S_c_aX)
        gV0_cS_aX.append(gr0_cS_aX)

        gr0_cI_aX=Jacob_IV(Qv,Qp, lat_per, ds,m,m2,  phiz, epso,beta, phi,V_aX[t], Z_aX[t], I_c_aX, S_aX[t])
        gV0_cI_aX.append(gr0_cI_aX)
        
        gr0_ncovIS_aX=Jacob_IV(Qv,Qp, lat_per, ds,m,m2, phiz,epso, beta, phi,V_aX[t], Z_aX[t], I_aX[t], S_aX[t_rand])
        gV0_ncovIS_aX.append(gr0_ncovIS_aX)
        
    epsZ_0_aX=np.mean(gZ0_aX)
    epsZ_S_aX=np.mean(gZ0_cI_aX)-epsZ_0_aX
    epsZ_I_aX=np.mean(gZ0_cS_aX)-epsZ_0_aX
    epsZ_IS_aX= np.mean(gZ_aX)-(epsZ_0_aX+epsZ_S_aX+epsZ_I_aX)
    epsZ_icovIS_aX=np.mean(gZ0_ncovIS_aX)-(epsZ_0_aX+epsZ_S_aX+epsZ_I_aX)
    epsZ_covIS_aX=epsZ_IS_aX-epsZ_icovIS_aX

    epsV_0_aX=np.mean(gV0_aX)
    epsV_S_aX=np.mean(gV0_cI_aX)-epsV_0_aX
    epsV_I_aX=np.mean(gV0_cS_aX)-epsV_0_aX
    epsV_IS_aX= np.mean(gV_aX)-(epsV_0_aX+epsV_S_aX+epsV_I_aX)
    epsV_icovIS_aX=np.mean(gV0_ncovIS_aX)-(epsV_0_aX+epsV_S_aX+epsV_I_aX)
    epsV_covIS_aX=epsV_IS_aX-epsV_icovIS_aX
    
    effs=[epsZ_0_aX, epsZ_S_aX, epsZ_I_aX, epsZ_icovIS_aX, epsZ_covIS_aX, epsV_0_aX, epsV_S_aX, 
          epsV_I_aX, epsV_icovIS_aX, epsV_covIS_aX]
    return effs
    
# Modern coexistence theory (simulations without Z or V and invasion analysis of the effects)
def MCT_analysis_SIVZ(mu, mui, lat_per, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, ntot):
    # absent/invading zoop
    init_conditions=[0.001, 0.001,0.001,0]
    result_az=simulation_SIVZ_rk4(mu, mui, lat_per, beta, phi, d, 
                               m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, init_conditions)
    
    Nt=len(result_az[0])
    i1=Nt-round(365*10/dt)
    i2=Nt
    S_az=result_az[0][i1:i2]
    I_az=result_az[1][i1:i2]
    V_az=result_az[2][i1:i2]
    Z_az=result_az[3][i1:i2]
    mZs_az=result_az[4][i1:i2]
    mVs_az=result_az[5][i1:i2]
    gZs_az=result_az[6][i1:i2]
    gVs_az=result_az[7][i1:i2]

    p1, p2= find_1_period(S_az)
    print('peaks az')
    print(p1)
    print(p2)
    S_az=S_az[p1:p2]
    I_az=I_az[p1:p2]
    Z_az=Z_az[p1:p2]
    V_az=V_az[p1:p2]

    S_c_az=np.mean(S_az)
    I_c_az=np.mean(I_az)
    
    surv_S_az=1
    if np.max(S_az)/Qp<1:
        surv_S_az=0
    if surv_S_az==1:
        effects_az=MCT_coexistence_effects_IS_on_VZ(S_az, I_az,Z_az,V_az, S_c_az, I_c_az, i1, i2,mu,eps,epso,d, phiz, eps_z, dz,dz2, m,m2, phi, Qp,Qv, beta, lat_per,alph,CC, 'zoop')
    else:
        effects_az=[np.nan for k in range(10)]
    # absent/invading virus
    init_conditions=[0.001, 0,0,0.001]
    result_av=simulation_SIVZ_rk4(mu, mui, lat_per, beta, phi, d, 
                               m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, init_conditions)

    S_av=result_av[0][i1:i2]
    I_av=result_av[1][i1:i2]
    V_av=result_av[2][i1:i2]
    Z_av=result_av[3][i1:i2]
    mZs_av=result_av[4][i1:i2]
    mVs_av=result_av[5][i1:i2]
    gZs_av=result_av[6][i1:i2]
    gVs_av=result_av[7][i1:i2]

    p1, p2= find_1_period(S_av)
    print('peaks av')
    print(p1)
    print(p2)
    S_av=S_av[p1:p2]
    I_av=I_av[p1:p2]
    Z_av=Z_av[p1:p2]
    V_av=V_av[p1:p2]

    S_c_av=np.mean(S_av)
    I_c_av=np.mean(I_av)
    
    surv_S_av=1
    if np.max(S_av)/Qp<1:
        surv_S_av=0
    if surv_S_av==1:
        effects_av=MCT_coexistence_effects_IS_on_VZ(S_av, I_av,Z_av,V_av, S_c_av, I_c_av, i1, i2,mu,eps,epso,d, phiz, eps_z, dz,dz2, m,m2, phi, Qp,Qv, beta, lat_per, alph, CC,'virus')
    else:
        effects_av=[np.nan for k in range(10)]

    delta_effects=[]
    for i in range(ntot-1):
        de=effects_av[i+ntot-1]-effects_av[i]
        delta_effects.append(de)
    for i in range(ntot-1):
        de0=effects_az[i]-effects_az[i+ntot-1]
        delta_effects.append(de0)
    return delta_effects

# find a period
def find_1_period(S):
    res=find_peaks(S)
    res=res[0]
    if len(res)>1:
        p1=res[0]
        p2=res[1]
        mx=max([S[p1], S[p2]])
        imx=np.argmax(np.array([S[p1], S[p2]]))
        inds=[p1,p2]
        if abs(S[p1]-S[p2]) > 0.1:
            smx=mx
            p_smx=inds[imx]
            if len(res)>2:
                for i in range(2, len(res)):
                    pn=res[i]
                    diff=abs(S[p_smx]-S[pn])
                    if diff<0.1:
                        p1=p_smx
                        p2=i
                        break
    else:
        # no period found => constant
        p1=0
        p2=1
    return p1, p2

# coexistence fitness analysis in the case of coexistence
def coexistence_analysis_SIVZ(S, I, V, Z,mu, mui, lp, beta, phi, d, m,m2, Qv, Qp,Qz,  eps, epso,phiz,eps_z,dz,dz2, CC,alph,dt, ndays, ntot):
    gZ=[]
    gZ0=[]
    gZ0_cS=[]
    gZ0_cI=[]
    gZ0_ncovIS=[]

    gV=[]
    gV0=[]
    gV0_cS=[]
    gV0_cI=[]
    gV0_ncovIS=[]

    nsteps=len(S)
    res=find_peaks(S)
    res=res[0]
    if len(res)>1:
        p1=res[0]
        p2=res[1]
        print('peaks found:')
        print(p1)
        print(p2)
        print(S[p1])
        print(S[p2])
        mx=max([S[p1], S[p2]])
        imx=np.argmax(np.array([S[p1], S[p2]]))
        inds=[p1,p2]
        if abs(S[p1]-S[p2]) > 0.1:
            smx=mx
            p_smx=inds[imx]
            if len(res)>2:
                for i in range(2, len(res)):
                    pn=res[i]
                    diff=abs(S[p_smx]-S[pn])
                    if diff<0.1:
                        p1=p_smx
                        p2=i
                        print('new peaks:')
                        print(p1)
                        print(p2)
                        print(S[p1])
                        print(S[p2])
                        break
    else:
        p1=0
        p2=len(S)
    S_c=np.mean(np.array(S[p1:p2]))
    I_c=np.mean(np.array(I[p1:p2]))
    for t in range(p1, p2):
        t_rand=random.randint(p1 , p2-1)
        gZ.append(phiz*eps_z*(S[t]+alph*I[t])-dz-dz2*Z[t])
        gZ0.append(phiz*eps_z*(S_c+alph*I_c)-dz-dz2*Z[t])
        gZ0_cS.append(phiz*eps_z*(S_c+alph*I[t])-dz-dz2*Z[t])
        gZ0_cI.append(phiz*eps_z*(S[t]+alph*I_c)-dz-dz2*Z[t])
        gZ0_ncovIS.append(phiz*eps_z*(S[t]+alph*I[t_rand])-dz-dz2*Z[t])
        
        
        gr=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S[t]/Qp-m2*V[t]
        gV.append(gr)
        
        gr0=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S_c/Qp-m2*V[t]
        gV0.append(gr0)
        
        gr0_cS=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S_c/Qp-m2*V[t]
        gV0_cS.append(gr0_cS)
        
        gr0_cI=beta/lp*(Qv/Qp)*I_c/V[t]-m-phi*S[t]/Qp-m2*V[t]
        gV0_cI.append(gr0_cI)
        
        gr0_ncovIS=beta/lp*(Qv/Qp)*I[t]/V[t]-m-phi*S[t_rand]/Qp-m2*V[t]
        gV0_ncovIS.append(gr0_ncovIS)

    epsZ_0=np.mean(gZ0)
    epsZ_S=np.mean(gZ0_cI)-epsZ_0
    epsZ_I=np.mean(gZ0_cS)-epsZ_0
    epsZ_IS= np.mean(gZ)-(epsZ_0+epsZ_S+epsZ_I)
    epsZ_icovIS=np.mean(gZ0_ncovIS)-(epsZ_0+epsZ_S+epsZ_I)
    epsZ_covIS=epsZ_IS-epsZ_icovIS

    epsV_0=np.mean(gV0)
    epsV_S=np.mean(gV0_cI)-epsV_0
    epsV_I=np.mean(gV0_cS)-epsV_0
    epsV_IS= np.mean(gV)-(epsV_0+epsV_S+epsV_I)
    epsV_icovIS=np.mean(gV0_ncovIS)-(epsV_0+epsV_S+epsV_I)
    epsV_covIS=epsV_IS-epsV_icovIS

    effs=[epsZ_0, epsZ_S, epsZ_I, epsZ_icovIS, epsZ_covIS, epsV_0, epsV_S,
          epsV_I, epsV_icovIS, epsV_covIS]
    
    return effs
