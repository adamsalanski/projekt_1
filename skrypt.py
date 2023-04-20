import numpy as np
from math import *

def Np(f,a,e2):
    N = a / np.sqrt(1-e2 * np.sin(f)**2) #**2 podnosi do kwadratu
    return(N)

def hirvonen(X,Y,Z,a,e2):
    p = np.sqrt(X**2 +Y**2)
    f = np.arctan(Z / (p * (1-e2))) # f to fi
    while True:
        N = Np(f,a,e2)
        h = (p/np.cos(f)) - N
        fp = f
        f = np.arctan(Z / (p * (1-e2 * N / (N +h))))
        if abs(fp - f) < (0.000001/206265):
            break
    l = np.arctan2(Y,X)    
    return(f,l,h)