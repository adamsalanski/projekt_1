import numpy as np
from math import *

class Transformacje:
    def __init__(self,model: str = "wgs84"):
        if model == "wgs84":
            self.a = 6378137.0
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "ewKrasowski":
            self.a = 6378245.0
            self.b = 6356863.019
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b)/self.a
        self.e2 = (2*self.flattening - self.flattening**2)
        print(model,self.a,self.b)
        
        def Rneu(self,f,l): # f to fi l to lambda
            R = np.array([[-np.sin(f) * np.cos(l),-np.sin(l),np.cos(f) * np.cos(l)],
                          [-np.sin(f) * np.sin(l),np.cos(l),np.cos(f) * np.sin(l)],
                          [np.cos(f), 0 ,np.sin(f)]])
            return(R)

        #zmiana XYZ na neu
        def XYZ2neu(self,dX,f,l):
            R = Rneu(f,l)
            return(R.T @ dX)
