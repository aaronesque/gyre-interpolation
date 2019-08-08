import numpy as np
import sympy as sym



def tau(xa,xb, x):
    h = xb-xa
    tau = (x-xa)/h
    return(tau)

def H30(t):
    H30 = 1 - 3*t**2 + 2*t**3
    return(H30)

def H31(t):
    H31 = t - 2*t**2 + t**3
    return(H31)

def H32(t):
    H32 = -t**2 + t**3
    return(H32)

def H33(t):
    H33 = 3*t**2 - 2*t**3
    return(H33)


def pressure(xa,xb, Pa,Pb, dPa,dPb):
    x = sym.Symbol('x')
    t = tau(xa,xb, x)
    h = xb-xa
    
    Pab = Pa*H30(t) + h*dPa*H31(t) + h*dPb*H32(t) + Pb*H33(t)
    
    return(Pab)


def dpressure(xa,xb, Pa,Pb, dPa,dPb):
    x = sym.Symbol('x')
    
    Pab = pressure(xa,xb, Pa,Pb, dPa,dPb)
    dPab = sym.diff( Pab, x )
    
    return(dPab)


def m2_integral(xa,xb, Pa,Pb, dPa,dPb):
    x = sym.Symbol('x')
    Pab = pressure(xa,xb, Pa,Pb, dPa,dPb)
    dPab = dpressure(xa,xb, Pa,Pb, dPa,dPb)
        
    m2ab = -8*np.pi*( sym.integrate( dPab*x**4 , (x, xa, xb) ) )

    return(m2ab)
