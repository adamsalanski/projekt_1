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
        


        def Np(self,f):
            N = self.a / np.sqrt(1-self.e2 * np.sin(f)**2) #**2 podnosi do kwadratu
            return(N)

        def hirvonen(self,X,Y,Z):
            p = np.sqrt(X**2 +Y**2)
            f = np.arctan(Z / (p * (1-self.e2))) # f to fi
            while True:
                N = Np(self,f)
                h = (p/np.cos(f)) - N
                fp = f
                f = np.arctan(Z / (p * (1-self.e2 * N / (N +h))))
                if abs(fp - f) < (0.000001/206265):
                    break
            l = np.arctan2(Y,X)    
            return(f,l,h)
        


   