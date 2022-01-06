#%%
"""
Created on Jan 2 2022
Main file: UT_sharedCode.py
Uncertain Finance
European options risks by uncertainty theory

@author: Carlos Alexander Grajales
Universidad de Antioquia, Medellín, Colombia
alexander.grajales@udea.edu.co
"""


"""
See working paper:  
Uncertainty and stochastic theories on European options valuation
and their delta and vega risks.
Grajales - Medina - Mongrut, January 1, 2022 
Available at SSRN: https://ssrn.com/abstract=3998282


Stochastic model, BSM: dX_t = r X_t dt + sigma X_t dW_t (1)
Uncertain model, Liu: dX_t = e X_t dt + sigma X_t dC_t  (2)
initial value X(0)=X_0

(a) Estimations of stock options: [sc,sp], [uc,up] 
     [sc,sp] stochastic call, stochastic put - from BSM model
     [uc,up] uncertain call, uncertain put - from Liu model
(b) Comparisons of stock options
     -norm metrics on matrices [uc_sc, up_sp]: Frobenius, norm1, norminf
(c) follow (a) and (b) for risk letters delta and vega
 
for paper we set:
we examine (a) and (b) at t = 1
X_0=1; 
grid sizes for e, sigma, K, T: nume=10; numsigma=10; numk=7; numt=4;
"""

import numpy as np
import scipy.stats as st
from scipy.integrate import quad_vec
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl
import enum

class OptionType(enum.Enum):
    CALL = 1.0
    PUT = -1.0

def BS_opt(S, K, r, q, sigma, T, CP, result):
    # result: 0->value; 1->delta; 2->gamma
       
    d       = (np.log(S / K) + \
                (r - q + 0.5 * np.power(sigma,2.0))* T) / \
                (1.0*sigma * np.sqrt(T))
    Nd1     = st.norm.cdf(d)
    Nd2     = st.norm.cdf(d - sigma*np.sqrt(T))
    BS_call = S * np.exp(-q*T)*Nd1 - K*np.exp(-r*T)*Nd2
    
    if CP == OptionType.CALL:
        price = BS_call
        delta = np.exp(-q*T)*Nd1
    elif CP == OptionType.PUT:
        price = BS_call - S*np.exp(-q*T) + K*np.exp(-r*T)
        delta = np.exp(-q*T)*(Nd1 - 1)
    
    Nd1prime = 1/np.sqrt(2*np.pi)*np.exp(-0.5*d*d)
    gamma    = Nd1prime*np.exp(-q*T) / (S*sigma*np.sqrt(T))
    
    if result == 0:
        BS_opt = price
    elif result == 1:
        BS_opt = delta
    elif result == 2:
        BS_opt = gamma
    return BS_opt

def Liu_opt(K,s,e,sigma,Y_0,r): #4D: K s e sigma
        Jc = lambda x: np.maximum( \
                    (Y_0*np.exp(e*s+(sigma*s*np.sqrt(3)/np.pi)* \
                    (np.log(x)-np.log(1-x))))-K,0)
        Jp = lambda x: np.maximum( \
                    K-(Y_0*np.exp(e*s+(sigma*s*np.sqrt(3)/np.pi)* \
                    (np.log(x)-np.log(1-x)))),0)
        minpy=1e-16
        uc, errc = quad_vec(Jc,minpy,1-minpy)
        up, errp = quad_vec(Jp,minpy,1-minpy)
        uc = np.exp(-r*s)*uc
        up = np.exp(-r*s)*up
        LiuOpt = {"uc":uc,"up":up }
        return LiuOpt

def mainCalculation():
    
    # Inputs
    
    X_0 = 1
    t   = 1
    Rt  = 0.01
    
    nume     = 10
    numsigma = 10
    v_e      = np.linspace(0.01,0.30,nume)
    v_sigma  = np.linspace(0.10,0.60,numsigma)
    sigma, e = np.meshgrid(v_sigma,v_e, indexing='xy')
    
    numk  = 7
    numt  = 4
    Str   = np.linspace(0.7*X_0,1.3*X_0,numk)
    Tm    = np.linspace(0.25*t,1*t,numt)
    
    # Stock options prices [sc,sp], [uc,up] 
    
        # stochastic stock options
        # we take variables: Strike, Time, Volatility (K,T,VOL)
    
    Vol = v_sigma
    sStr,sVol,sTm = np.meshgrid(Str, Vol, Tm) #3D arrange as (K x T) x VOL !!!
    
    CP  = OptionType.CALL
    sc  = BS_opt(X_0,sStr,Rt,0,sVol, sTm, CP, 0)
    CP  = OptionType.PUT
    sp  = BS_opt(X_0,sStr,Rt,0,sVol, sTm, CP, 0)       
    
        # uncertain stocks options
        # we take variables: Strike, Time, Mean, Volatility (K,T,Mu,VOL)
    
    Mn = v_e;
    [uTm,uStr,uMn,uVol] = np.meshgrid(Tm,Str,Mn,Vol);
    
    LiuOpt = Liu_opt(uStr,uTm,uMn,uVol,X_0,Rt)
    uc = LiuOpt["uc"]
    up = LiuOpt["up"]
    
    #Figure
    #rstride=40, cstride=40
    kk, tt = np.meshgrid(Str,Tm, indexing='ij')

    ax = plt.axes(projection='3d')
    ax.plot_wireframe(kk, tt, sc[5,:,:], rstride=40, cstride=40)
    ax.plot_surface(kk, tt, uc[:,:,1,1], color='m', alpha = .5)
    ax.plot_surface(kk, tt, uc[:,:,5,5], color='r', alpha = .5)
    ax.set_yticks([0.25, 0.50, 0.75, 1.00])
    ax.set_xlabel("strike")
    ax.set_ylabel("maturity")
    ax.set_zlabel("call price")
    ax.set_title("BSM (stoch.) vs LIU (uncert.)")
    bsm = mpl.lines.Line2D([0],[0], linestyle="none", c='b', marker = 'o')
    liu1 = mpl.lines.Line2D([0],[0], linestyle="none", c='m', marker = 'o')
    liu2 = mpl.lines.Line2D([0],[0], linestyle="none", c='r', marker = 'o')
    ax.legend([bsm, liu1, liu2], ['BSM', 'Liu 1', 'Liu 2'], numpoints = 1, loc='center left', bbox_to_anchor=(1.17, 0.5), fontsize=7)
    

    print("BSM surface, sigma = {0:10.5f}".format(Vol[5]))
    print("Liu surface 1, expert e = {0:10.5f}, expert sigma = {1:10.5f}".format(Mn[1], Vol[1]))
    print("Liu surface 2, expert e = {0:10.5f}, expert sigma = {1:10.5f}".format(Mn[5], Vol[5]))
    
    return {"sc":sc,"sp":sp,"uc":uc,"up":up}

    
y=mainCalculation()