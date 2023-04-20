import numpy as np
from math import *

def Rneu(f,l): # f to fi l to lambda
    R = np.array([[-np.sin(f) * np.cos(l),-np.sin(l),np.cos(f) * np.cos(l)],
                  [-np.sin(f) * np.sin(l),np.cos(l),np.cos(f) * np.sin(l)],
                  [np.cos(f), 0 ,np.sin(f)]])
    return(R)

#zmiana XYZ na neu
def XYZ2neu(dX,f,l):
    R = Rneu(f,l)
    return(R.T @ dX)