
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
        
        def BL22000 (self, f, l):
            if (l >13.5 and  l < 16.5):
                zone =5
                l0 = 15
            elif (l > 16.5 and l < 19.5):
                zone = 6 
                l0 = 18
            elif (l > 19.5 and l < 22.5):
                zone =7
                l0 = 21
            elif (l > 22.5 and l <25.5):
                zone = 8
                l0 = 24
            f = f * pi/180
            l = l * pi/180
            l0 = l0 * pi/180
            b2 = (self.a**2) * (1-self.e2)
            ep2 = (self.a**2-b2)/(b2)
            t = atan(f)
            n2 = ep2 *((cos(f))**2)
            N = self.a / (sqrt(1 - self.e2 * (sin(f)) ** 2))
            A0 = 1-(self.e2/4) - ((3*(self.e2**2))/64)-((5*(self.e2**3))/256)
            A2 = (3/8)*(self.e2+(self.e2**2/4)+((15*(self.e2**3))/128))
            A4 = (15/256)*((self.e2**2)+((3*(self.e2**3))/4))
            A6 = (35*(self.e2**3))/3072
            sigma = self.a * (A0*(f) - A2 * sin(2*f) + A4 * sin(4*f) - A6 * sin(6 * f))
            dlam = l - l0
            xgk = sigma + ((dlam**2)/2) * N * sin(f) * cos(f) * (1 + ((dlam**2)/12) * ((cos(f))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((dlam**4)/360) * ((cos(f))**4) * (61-58 * (t**2) + t**4 + 270 * n2 - 330 * n2 * (t**2)))
            ygk = dlam * N * cos(f) * (1 + ((dlam**2)/6) * ((cos(f))**2) * (1 - t**2 + n2) + ((dlam**4)/120) * ((cos(f))**4) * (5 - 18 * (t**2) + t**4 + 14 * n2 - 58 * n2 * (t**2)))
            m = 0.999923 #skala PL-2000
            x2000 = xgk * m
            y2000 = ygk * m + (zone * 1000000) + 500000
            return x2000, y2000

            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
