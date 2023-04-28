import numpy as np
from math import *
import argparse




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
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie dziesiętne]

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
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie dziesiętne]
        l : TYPE: [float] - Długoć geodezyjna [stopnie dziesiętne]
            DESCRIPTION.

        Returns
        -------
        x92 : TYPE: [float] - Współrzędna X w układzie PL-1992 [metry]
        y92 : TYPE: [float] - Współrzędna Y w układzie PL-1992 [metry]

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
        f : TYPE: [float] - Szerokoć geodezyjna [stopnie dziesiętne]

        Returns
        -------
        N : TYPE: [float] - Promień krzywizny w I wertykale [metry]

        '''
        N = self.a / np.sqrt(1-self.e2 * np.sin(f)**2) #**2 podnosi do kwadratu
        return(N)
    
    def BL22000 (self, f, l,ns):
        '''
        Funkcja przeliczająca współrzędne geodezyjne na współrzędne w układzie 2000.
        
        Parameters
        ----------
        f : TYPE : [float] : Szerokość geodezyjna [stopnie dziesiętne]
        l : TYPE : [float] : Długość geodezyjna [stopnie dziesiętne]
        ns : TYPE: [int] : Numer strefy

        Returns
        -------
        x2000 : TYPE : [float] : współrzędna X w układzie 2000 [metry]
        y2000 : TYPE : [float] : współrzędna Y w układzie 2000 [metry]

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
    
    def Rneu(self, f, l):
        """
        Funkcja, która, przyjmujac współrzedne krzywoliniowe utworzy macierz obrotu 
        potrzebną do przeliczenia współrzędnych do układu współrzędnych neu
    
        INPUT:
        ----------
        f : [float] : wspołrzędna fi punktu początkowego układu lokalnego [stopnie dziesiętne]
        l : [float] :wspołrzędna l punktu początkowego układu lokalnego [stopnie dziesiętne]
    
        OUTPUT:
        -------
        R : [array of float64] : macierz obrotu
    
        """
        N=[(-sin(f) * cos(l)), (-sin(f) * sin(l)), (cos(f))]
        E=[(-sin(l)), (cos(l)),  (0)]
        U=[( cos(f) * cos(l)), ( cos(f) * sin(l)), (sin(f))]
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
    
    
    def XYZ2neu(self,s,alfa,z,X,Y,Z):
        '''
        Funkcja przeliczająca współrzędzne kartezjańskie X,Y,Z do układu NEU przy pomocy wektora przestrzennego
        w układzie topocentrycznym.

        Parameters
        ----------
        Parametry wektora przestrzennego w układzie lokalnym topocentrycznym:
        s : TYPEE: [float] - Odległosć [metry]
        alfa : TYPE: [float] - Kąt pionowy [stopnie dziesiętne]
        z : TYPE:[float] - Kąt poziomy [stopnie dziesiętne]
        X : TYPE: [float] - Współrzędna X w układzie ortokartezjańskim [metry]
        Y : TYPE: [float] - Współrzędna Y w układzie ortokartezjańskim [metry]
        Z : TYPE: [float] - Współrzędna Z w układzie ortokartezjańskim [metry]

        Returns
        -------
        neu : TYPE: [array of float64] - współrzedne topocentryczne (North , East (E), Up (U)) 

        '''
        dX = np.array([s * np.sin(z) * np.cos(alfa),# to a to alfa
                         s * np.sin(z) * np.sin(alfa),# to a to alfa
                         s * np.cos(z)])
        f,l,h = self.hirvonen(X,Y,Z)
        
        N=[(-sin(f) * cos(l)), (-sin(f) * sin(l)), (cos(f))]
        E=[(-sin(l)), (cos(l)),  (0)]
        U=[( cos(f) * cos(l)), ( cos(f) * sin(l)), (sin(f))]
        R=np.transpose(np.array([N,E,U]))
        neu = R.T @ dX
        return(neu)
    
    def zargparse(self):
            
            parser = argparse.ArgumentParser(description='Transformacje wspolrzednych')
            
            parser.add_argument('-X', help='wartosc wspolrzednej pierwszegi punktu X [m]', required=False, type=float)
            parser.add_argument('-Y', help='wartosc wspolrzednej pierwszego punktu Y [m]', required=False, type=float)
            parser.add_argument('-Z', help='wartosc wspolrzednej pierwszego punktu Z [m]', required=False, type=float)
               
            parser.add_argument('-f', help="wartosc wspolrzednej f [Â° ' '']", required=False, type=float)
            parser.add_argument('-l', help="wartosc wspolrzednej l [Â° ' '']", required=False, type=float)
            parser.add_argument('-H', help='wartosc wspolrzednej H [m]', required=False, type=float)
            
            parser.add_argument('-s', help='wartosc dlugosci miedzy dwoma punktami [m]', required=False, type=float)
            parser.add_argument('-alfa', help="wartosc kat poziomego [Â° ' '']", required=False, type=float)
            parser.add_argument('-z', help="wartosc kat zenitalnego [Â° ' '']", required=False, type=float)
            parser.add_argument('-ns', help="numer strefy w odwzorowaniu PL-2000 [Â° ' '']", required=False, type=float)
            parser.add_argument('-model', help='model elipsoidy', choices=['grs80','wgs84', 'ewKrasowski'], required=False, type=str, default='grs80')
            args = parser.parse_args()
            
            self.__init__(args.model)
            
            return( args.X, args.Y, args.Z, args.f, args.l, args.H, args.s, args.alfa, args.z, args.ns)

if __name__ == "__main__":
    geo = Transformacje(model = "grs80")
    
    print('Przykładowe wywołanie hirvonen')
    X = 3664940.500
    Y = 1409153.590
    Z = 5009571.170
    phi, lam, h = geo.hirvonen(X, Y, Z)
    print(phi, lam, h)
    
    print('Przykładowe wywołanie flh2XYZ')
    X,Y,Z = geo.flh2XYZ(phi,lam,h)
    print(X,Y,Z)
    
    print('Przykładowe wywołanie fl2pl1992')
    f=0.9076010398716878
    l = 0.27936573137624443 
    x92,y92 = geo.fl2pl1992(f, l)
    print(x92,y92)
    
    print('Przykładowe wywołanie BL22000')
    x2000,y2000 = geo.BL22000(f,l,5)
    print(x2000,y2000)
    
    print('Przykładowe wywołanie XYZ2neu')
    s = 30000.000 
    alfa = 165.0000 
    z = 90.0000
    neu = geo.XYZ2neu(s, alfa, z, X, Y, Z)
    print(neu)
    
    plik = "wsp_inp.txt"
    tablica = np.genfromtxt(plik, delimiter=',', skip_header = 4)
    dane = np.ones((12,3))
    for i, n in enumerate(tablica):
        phi, lam, hel = geo.hirvonen(tablica[i,0], tablica[i,1], tablica[i,2])
        dane[i,:] = [phi, lam, hel]
    print('f,l,h zapis do pliku txt')    
    print(dane) #phi, lam, hel
    
    dane2 = np.ones((12,3))
    for i, n in enumerate(dane):
        X, Y, Z = geo.flh2XYZ(dane[i,0], dane[i,1], dane[i,2])
        dane2[i,:] = [X, Y, Z]
    print('X,Y,Z zapis do pliku txt')     
    print(dane2) # X, Y, Z
    
    
        
    dane3 = np.ones((12,2))
    for i, n in enumerate(dane):
        x2000, y2000 = geo.BL22000(dane[i,0], dane[i,1],7)
        dane3[i,:] = [x2000, y2000]
    print('x2000,y2000 zapis do pliku txt')     
    print(dane3) # x2000, y2000
    
    dane4 = np.ones((12,2))
    for i, n in enumerate(dane):
        x1992, y1992 = geo.fl2pl1992(dane[i,0], dane[i,1])
        dane4[i,:] = [x1992, y1992]
    print('x1992,y1992 zapis do pliku txt')     
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
    print('neu zapis do pliku txt') 
    print(neu) 
    
    DANE = np.hstack((dane, dane2, neu, dane3, dane4))
    np.savetxt("wyniki.txt", DANE, delimiter='  ', fmt = ['%10.8f', '%10.8f', '%10.5f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f', '%10.8f' ], header = 'Konwersja współrzednych geodezyjnych ', comments=' phi [st]     lambda[st]     hel[m]          X[m]              Y[m]              Z[m]          N[m]         E[m]         U[m]         x2000[m]        y2000[m]        x1992[m]          y1992[m]      \n ' )
    
    

    proba = Transformacje()
    X,Y,Z,f,l,h,s,alfa,z,ns = proba.zargparse()
    print('_______________________________________________')
    print('___________DANE:_______________________________')
    print('(X,Y,Z)',X,Y,Z)
    print('(f,l,h)',f,l,h)
    print('(s,alfa,z)',s, alfa, z)
    print('(ns)',ns)
    print('_______________________________________________')
    print('__________WYNIKI:______________________________')
    try:
        f1,l1,h1 = proba.hirvonen(X,Y,Z) 
        print('(f,l,h)',f1,l1,h1)
    except TypeError: 
        pass
    try:
        X2,Y2,Z2 = proba.flh2XYZ(f,l,h)
        print('(X,Y,Z)',X2,Y2,Z2)
    except TypeError:
        pass
    try:
        x92,y92 = proba.fl2pl1992(f, l)
        print('(x92,y92)',x92,y92)
    except TypeError:
        pass
    try:
        x2000,y2000 = proba.BL22000(f, l, ns)
        print('(x2000,y2000)',x2000,y2000)
    except TypeError:
        pass
    except UnboundLocalError:
        pass
    try:
        neu = proba.XYZ2neu(s, alfa, z, X, Y, Z)
        print('neu',neu)
    except TypeError:
        pass
    
    
    
        

    
    
            
            
            
    

    
    
    
    

    
    
     
     
     
     
     
    


        

        


   

        
        
         
         



