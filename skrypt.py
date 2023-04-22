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
        
    def sigma(self,f):
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigma = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        return sigma
    
    def fl2pl1992(self, f,l):
        f = f * pi/180
        l = l * pi/180
        m0 = 0.9993
        l0 = 19 * pi/180
        b2 = self.a**2*(1 - self.e2)
        ep2 = (self.a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = Np(self,f)
        sigm = sigma(self,f)
        xgk = sigm + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x92 = xgk * m0 - 5300000
        y92 = ygk * m0 + 500000
        return x92,y92


        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            

