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
    
    def dms(self,x):
        '''
        Funkcja zamieniająca stopnie dziesiętne na stopnie, minuty, sekundy.

        Parameters
        ----------
        x : TYPE: [float] - Kąt [stopnie dziesiętne]

        Returns
        -------
        d,m,s - Kąt [stopnie]
        
        d : TYPE: [int] -  [stopnie]
        m : TYPE: [int] -  [minuty]
        s : TYPE: [float] - [sekundy]

        '''
        znak = ' '
        if x<0:
            znak = '-'
            x = abs(x)
        x = x * 180/pi
        d = int(x) #stopnie
        m = int((x - d)*60) #minuty
        s = (x - d - m/60)*3600 #sekundy
        print(znak,"%3d %2d %7.5f"% (d,m,s)) 
        return(d,m,s)
    
    def fromdms(X):
        '''
        Funkcja zamieniająca stopnie w układzie (d m s) na radiany i stopnie dziesiętne.
         
        Parameters
        ----------
        X : Kąt wyrażony w postaci d m s, gdzie
        d : TYPE: [int] -  [stopnie]
        m : TYPE: [int] -  [minuty]
        s : TYPE: [float] - [sekundy]

        Returns
        -------
        Z : TYPE: [float] - Kąt [radiany]
        Y : TYPE: [float] - Kąt [stopnie dziesiętne]

        '''
      
        znak = 1
        if X[0] == '-':
             znak = -1
        Y = X.split(' ' or '"' or "'" or '°' or '-')
        d = int(Y[0])
        m = int(Y[1])
        s = float(Y[2])
        s = s/3600
        m = m/60
        Y = znak*(d+m+s)
        Z = Y * pi/180
        Y = float(f'{Y:7.5f}')
        return(Z,Y)
        


    def sigma(self,f):
        '''
        Funkcja licząca sigmę wykorzystywaną w funkcjach fl2pl1992 oraz BL2200.

        Parameters
        ----------
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie]

        Returns
        -------
        sigma : TYPE: [float] - Sigma

        '''
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigma = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        return sigma
    
    def fl2pl1992(self, f,l):
        '''
        Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 1992.

        Parameters
        ----------
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie]
        l : TYPE: [float] - Długoć geodezyjna [stopnie]
            DESCRIPTION.

        Returns
        -------
        x92 : TYPE: [float] - Współrzędna X w układzie PL-1992 [m]
        y92 : TYPE: [float] - Współrzędna Y w układzie PL-1992 [m]

        '''
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

    def Np(self,f):
        '''
        Funkcja licząca promień krzywizny w I wertykale.

        Parameters
        ----------
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie]

        Returns
        -------
        N : TYPE: [float] - Promień krzywizny w I wertykale [m]

        '''
        N = self.a / np.sqrt(1-self.e2 * np.sin(f)**2) #**2 podnosi do kwadratu
        return(N)
    
    def BL22000 (self, f, l):
        '''
        Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 2000.
        
        Parameters
        ----------
        f : [float] : Szeroko
        l : TYPE
            DESCRIPTION.

        Returns
        -------
        x2000 : TYPE
            DESCRIPTION.
        y2000 : TYPE
            DESCRIPTION.

        '''
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
    
    def hirvonen(self,X,Y,Z):
        '''
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.

        Parameters
        ----------
        X : TYPE: [float] - Współrzędna X w układzie ortokartezjańskim [metry]
        Y : TYPE: [float] - Współrzędna Y w układzie ortokartezjańskim [metry]
        Z : TYPE: [float] - Współrzędna Z w układzie ortokartezjańskim [metry]

        Returns
        -------
        f : TYPE: [float] - Szerokość geodezyjna [stopnie dziesiętne]
        l : TYPE: [float] - Długoć geodezyjna [stopnie dziesiętne]
        h : TYPE: [float] - Wysokoć elipsoidalna [metry]

        '''
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
    
    def flh2XYZ(f,l,h,self):
        '''
        Funkcja przelicza współrzędne krzywoliniowe(f,l,h) na współrzędne prostokątne(X,Y,Z).

        Parameters
        ----------
        f : TYPE: [float] - Szerokość geodezyjna [stopnie dziesiętne]
        l : TYPE: [float] - Długosć geodezyjna [stopnie dziesiętne]
        h : TYPE: [float] - Wysokosc elipsoidalna [metry]

        Returns
        ----------
        X : TYPE: [float] - Współrzędna X w układzie ortokartezjańskim [metry]
        Y : TYPE: [float] - Współrzędna Y w układzie ortokartezjańskim [metry]
        Z : TYPE: [float] - Współrzędna Z w układzie ortokartezjańskim [metry]
        

        '''
        f = f * pi/180
        l = l * pi/180
        N = Np(self,f)
        X = (N + h) * cos(f) * cos(l) 
        Y = (N + h) * cos(f) * sin(l)
        Z = (N + h - N * self.e2) * sin(f)
        return(X,Y,Z)
     
    def Rneu(self,f,l):
        '''
        Funkcja przyjmuje współrzedne krzywoliniowe i tworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych krzywoliniowych (f,l) do układu współrzędnych neu

        Parameters
        ----------
        f : TYPE: [float] - Szerokość geodezyjna [stopnie dziesiętne]
        l : TYPE: [float] - Długosć geodezyjna [stopnie dziesiętne]

        Returns
        -------
        R : [array of float64] : Macierz obrotu

        '''
        R = np.array([[-np.sin(f) * np.cos(l),-np.sin(l),np.cos(f) * np.cos(l)],
                       [-np.sin(f) * np.sin(l),np.cos(l),np.cos(f) * np.sin(l)],
                       [np.cos(f), 0 ,np.sin(f)]])
        return(R)

     
    def XYZ2neu(self,dX,f,l):
        '''
        Funckja obliczająca wektor w układzie neu

        Parameters
        ----------
        dX : TYP.
        f : TYPE: [float] - Szerokość geodezyjna [stopnie dziesiętne]
        l : TYPE: [float] - Długosć geodezyjna [stopnie dziesiętne]

        Returns
        
        R.T @ dX : [array of float64] : współrzedne topocentryczne (North , East (E), Up (U))

        '''
        R = Rneu(f,l)
        return(R.T @ dX)
     
    
     
     
     
     
     
    


        

        


   

        
        
         
         



