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

