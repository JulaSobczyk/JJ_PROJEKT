from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np

o = object()

class Transformacje:
    def __init__(self, model: str = "wgs84"):
        """
        Parametry elipsoid:
            a - duża półoś elipsoidy - promień równikowy
            b - mała półoś elipsoidy - promień południkowy
            flat - spłaszczenie
            ecc2 - mimośród^2
        + WGS84: https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84
        + Inne powierzchnie odniesienia: https://en.wikibooks.org/wiki/PROJ.4#Spheroid
        + Parametry planet: https://nssdc.gsfc.nasa.gov/planetary/factsheet/index.html
        """
        if model == "wgs84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "mars":
            self.a = 3396900.0
            self.b = 3376097.80585952
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
    # def Np(self, fi):     #promien krzywizny w I wertykale
    #     N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
    #     return(N)
    
        
    def xyz2plh(self, X, Y, Z, jedn = 'dec_degree'):
        """
        Algorytm Hirvonena - algorytm transformacji współrzędnych ortokartezjańskich (x, y, z)
        na współrzędne geodezyjne długość szerokość i wysokośc elipsoidalna (phi, lam, h). Jest to proces iteracyjny. 
        W wyniku 3-4-krotneej iteracji wyznaczenia wsp. phi można przeliczyć współrzędne z dokładnoscią ok 1 cm.     
        Parameters
        ----------
        X, Y, Z : FLOAT
             współrzędne w układzie orto-kartezjańskim, 

        Returns
        -------
        lat
            [stopnie dziesiętne] - szerokość geodezyjna
        lon
            [stopnie dziesiętne] - długośc geodezyjna.
        h : TYPE
            [metry] - wysokość elipsoidalna
        output [STR] - optional, default 
            dec_degree - decimal degree
            dms - degree, minutes, sec
        """
        r   = sqrt(X**2 + Y**2)           # promień
        phi_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        phi = 0
        while abs(phi_prev - phi) > 0.000001/206265:    
            phi_prev = phi
            N = self.a / sqrt(1 - self.ecc2 * sin(phi_prev)**2)
            h = r / cos(phi_prev) - N
            phi = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lam = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(phi))**2);
        h = r / cos(phi) - N       
        if jedn == "dec_degree":
            return degrees(phi), degrees(lam), h 
        elif jedn == "dms":
            phi = self.deg2dms(degrees(phi))
            lam = self.deg2dms(degrees(lam))
            return f"{phi[0]:02d}:{[1]:02d}:{phi[2]:.2f}", f"{lam[0]:02d}:{lam[1]:02d}:{lam[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{jedn} - jednostka niezdefiniowana")
        
    
    def plh2xyz(self, phi, lam, h):
        phi = radians(phi)
        lam = radians(lam)
        Rn = self.a / sqrt(1 - self.ecc2 * sin(phi)**2)
        q = Rn * self.ecc2 * sin(phi)
        x = (Rn + h) * cos(phi) * cos(lam)
        y = (Rn + h) * cos(phi) * sin(lam)
        z = (Rn + h) * sin(phi) * q
        return x, y, z
    
    