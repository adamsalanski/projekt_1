import numpy as np
from math import *

def Np(f,a,e2):
    N = a / np.sqrt(1-e2 * np.sin(f)**2) #**2 podnosi do kwadratu
    return(N)

def flh2XYZ(f,l,h,a,e2):
     N = Np(f,a,e2)
     X = (N + h) * cos(f) * cos(l) 
     Y = (N + h) * cos(f) * sin(l)
     Z = (N + h - N * e2) * sin(f)
     return(X,Y,Z)
 

