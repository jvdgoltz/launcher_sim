import math
import isatmos as isa
Ra = 8314.

def Tcb(T_f0,Mw,kappa,OFR,Hv):
    '''
    This is not a correct estimation based on chemical thermodynamics!
    '''
    Tc = T_f0 + Hv/((kappa*Ra)/(Mw*(kappa-1)) *(1+OFR))
    return Tc

def Tcm(T_f0,Mw,kappa,Hd):
    '''
    This is not a correct estimation based on chemical thermodynamics!
    '''
    Tc = T_f0 + Hd / ((kappa*Ra)/(Mw*(kappa-1)))
    return Tc

def vdk(kappa):
    '''
    returns Vandenkerckhove function
    for propellant specific heat ratio kappa
    '''
    return math.sqrt(kappa*((kappa+1)/2)**((1+kappa)/(1-kappa)))

    
def mdot(At,pc,Tc,Mw,kappa):
    '''
    returns engine mass flow;
    Input:
    throat area At,
    chamber pressure pc,
    chamber temperature Tc,
    molecular weight Mw,
    propellant specific heat ratio kappa
    '''    
    return pc*At/(math.sqrt((Ra/Mw)*Tc)) * vdk(kappa)
    
def vexit(epsilon,Tc,Mw,kappa):
    '''
    returns engine exit velocity;
    Input:
    area expansion ratio epsilon,
    chamber temperature Tc,
    molecular weight Mw,
    propellant specific heat ratio kappa
    '''
    pr = pratio(epsilon,kappa)
    vexit = math.sqrt(2*kappa/(kappa-1) * Ra/Mw * Tc *(1-pr**((kappa-1)/kappa)) )
    return vexit
    
def pratio(epsilon,kappa):
    '''
    returns ratio between nozzle exit pressure and burning chamber pressure;
    Input:
    area expansion ratio epsilon,
    propellant specific heat ratio kappa,
    Newton's Method is used
    '''
    Gamma = vdk(kappa)
    pratio = 0                                                                  #initial guess
    f = (Gamma**2/((2*kappa/(kappa-1))*epsilon**2))**(kappa/2)
    while f<>0:
        f = (Gamma**2/((2*kappa/(kappa-1))*epsilon**2) + pratio**((1+kappa)/kappa))**(kappa/2) - pratio                                                                 #calculate f(x)
        fp = ((1+kappa)/kappa) * pratio**((1+kappa)/kappa -1) * (kappa/2) * (Gamma**2/((2*kappa/(kappa-1))*epsilon**2) + pratio**((1+kappa)/kappa))**(kappa/2-1) - 1    #calculate f'(x)
        pratio=pratio - f/fp                                                    #Newton's Method x_(n+1) = x_n - f(x_n) / f'(x_n)
    return pratio
    
def FT(mdot,ve,At,epsilon,pc,pratio,pa):
    '''
    returns engine thrust:
    '''
    FT = mdot*ve + At*epsilon*(pc*pratio-pa)
    return FT