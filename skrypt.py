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
    
    def fromdms(self,X):
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
        
        m0 = 0.9993
        l0 = 19 * pi/180
        b2 = self.a**2*(1 - self.e2)
        ep2 = (self.a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.a / (sqrt(1 - self.e2 * (sin(f)) ** 2))
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigma = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sigma + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
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
    
    def BL22000 (self, f, l,ns):
        '''
        Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 2000.
        
        Parameters
        ----------
        f : TYPE : [float] : Szerokość geodezyjna [stopnie]
        l : TYPE : [float] : Długość geodezyjna [stopnie]
        ns : TYPE: [int] : Numer strefy

        Returns
        -------
        x2000 : TYPE : [float] : współrzędna X w układzie 2000 [m]
        y2000 : TYPE : [float] : współrzędna Y w układzie 2000 [m]

        '''
        if ns == 5:
            l0 = radians(15)
        elif ns == 6:
            l0 = radians(18)
        elif ns == 7:
            l0 = radians(21)
        elif ns == 8:
            l0 = radians(24)
        m0= 0.999923
        b2 = self.a**2*(1 - self.e2)
        ep2 = (self.a**2 - b2)/b2
        dl = l - l0
        t = tan(f)
        n2 = ep2 * cos(f)**2
        N = self.a / (sqrt(1 - self.e2 * (sin(f)) ** 2))
        A0 = 1 - self.e2/4 - 3 * self.e2**2/64 - 5 * self.e2**3/256
        A2 = (3/8) * (self.e2 + self.e2**2/4 + 15*self.e2**3/128)
        A4 = (15/256) * (self.e2**2 + (3 * self.e2**3)/4)
        A6 = 35 * self.e2**3/3072
        sigma = self.a * (A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))
        xgk = sigma + (dl**2/2) * N * sin(f)*cos(f)*(1 + (dl**2/12)*cos(f)**2*(5-t**2+9*n2+4*n2**2)+ ((dl**4)/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dl*N*cos(f)*(1+(dl**2/6)*cos(f)**2*(1 - t**2 + n2) + (dl**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x2000 = xgk * m0
        y2000 = ygk * m0 + ns * 1000000 + 500000
        return x2000,y2000
    
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
            N = self.a / np.sqrt(1-self.e2 * np.sin(f)**2)
            h = (p/np.cos(f)) - N
            fp = f
            f = np.arctan(Z / (p * (1-self.e2 * N / (N +h))))
            if abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)    
        return(f,l,h)
    
    def flh2XYZ(self,f,l,h):
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
        
        N = self.a / np.sqrt(1-self.e2 * np.sin(f)**2)
        X = (N + h) * cos(f) * cos(l) 
        Y = (N + h) * cos(f) * sin(l)
        Z = (N + h - N * self.e2) * sin(f)
        return(X,Y,Z)
    
    def średnia(self, wartosci):
        """
        Funkcja liczy średnią wartość z elementów w liscie
        
        Parameters:
        ----------
        wartosci : [float] : lista wartosci
        
        Returns:
        --------
        srednia : [float] : średnia arytmetyczna elementów z listy 
        
        """
        suma = 0
        ilość = 0
        for wartosc in wartosci:
            suma += wartosc
            ilość += 1
        srednia = float(suma / ilość)
        return(srednia)
    
    def Rneu(self, phi, lam):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        phi : [float] : wspołrzędna fi punktu początkowego układu lokalnego
        lam : [float] :wspołrzędna l punktu początkowego układu lokalnego
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        N=[(-sin(phi) * cos(lam)), (-sin(phi) * sin(lam)), (cos(phi))]
        E=[(-sin(lam)), (cos(lam)),  (0)]
        U=[( cos(phi) * cos(lam)), ( cos(phi) * sin(lam)), (sin(phi))]
        R=np.transpose(np.array([N,E,U]))
        return (R, N, E, U)
    
    def NEU(self, R,v):
        """
        Funckja obliczająca wektor w układzie neu
    
        Parameters:
        -----------
        R : R : [array of float64] : macierz obrotu
        v : [array of float64] : wektor w układzie XYZ
        
        Returns:
        -------
        NEU : [array of float64] : współrzedne topocentryczne (North , East (E), Up (U))
    
        """
        NEU=np.zeros(v.shape)
        for a in range(v.shape[0]):
            for b in range(3):
                for c in range(3):
                    NEU[a,c]+=v[a,b]*R[c,b]
        return (NEU)


if __name__ == "__main__":
    geo = Transformacje(model = "grs80")
    X = 3664940.500
    Y = 1409153.590
    Z = 5009571.170
    phi, lam, h = geo.hirvonen(X, Y, Z)
    print(phi, lam, h)
    
    X,Y,Z = geo.flh2XYZ(phi,lam,h)
    print(X,Y,Z)
    
    f=0.9076010398716878
    l = 0.27936573137624443 
    x92,y92 = geo.fl2pl1992(f, l)
    print(x92,y92)
    
    x2000,y2000 = geo.BL22000(f,l,5)
    print(x2000,y2000)
    
    plik = "wsp_inp.txt"
    tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
    dane = np.ones((12,3))
    for i, n in enumerate(tablica):
        phi, lam, hel = geo.hirvonen(tablica[i,0], tablica[i,1], tablica[i,2])
        dane[i,:] = [phi, lam, hel]
        
    print(dane) #phi, lam, hel
    
    dane2 = np.ones((12,3))
    for i, n in enumerate(dane):
        X, Y, Z = geo.flh2XYZ(dane[i,0], dane[i,1], dane[i,2])
        dane2[i,:] = [X, Y, Z]
        
    print(dane2) # X, Y, Z
    
    
        
    dane3 = np.ones((12,2))
    for i, n in enumerate(dane):
        x2000, y2000 = geo.BL22000(dane[i,0], dane[i,1],7)
        dane3[i,:] = [x2000, y2000]
        
    print(dane3) # x2000, y2000
    
    dane4 = np.ones((12,2))
    for i, n in enumerate(dane):
        x1992, y1992 = geo.fl2pl1992(dane[i,0], dane[i,1])
        dane4[i,:] = [x1992, y1992]
        
    print(dane4) # x1992, y1992
    
    phi_sr = geo.średnia(dane[:,0])
    lam_sr = geo.średnia(dane[:,1]) 
    X_sr = geo.średnia(dane2[:,0])
    Y_sr = geo.średnia(dane2[:,1])
    Z_sr = geo.średnia(dane2[:,2])
    [phi_sr, lam_sr, hel_sr] = geo.hirvonen(X_sr, Y_sr, Z_sr)
    R, N, E, U = geo.Rneu(phi_sr,lam_sr)
    
    ii=0
    v=np.array(np.zeros((dane2.shape[0],3)))
    for ii in range(0, dane2.shape[0]):
        v[ii,0]=X_sr - dane2[ii,0]
        v[ii,1]=Y_sr - dane2[ii,1]
        v[ii,2]=Z_sr - dane2[ii,2]
        ii += 1

    neu=geo.NEU(R, v)
    print(neu) 
    
    DANE = np.hstack((dane, dane2, neu, dane3, dane4))
    np.savetxt("wyniki.txt", DANE, delimiter='  ', fmt = ['%10.8f', '%10.8f', '%10.5f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f' ], header = 'Konwersja współrzednych geodezyjnych ', comments=' phi [st]     lambda[st]     hel[m]          X[m]              Y[m]              Z[m]          N[m]         E[m]         U[m]         x2000[m]        y2000[m]        x1992[m]          y1992[m]      \n ' )
    
    

    
    
     
     
     
     
     
    


        

        


   

        
        
         
         



