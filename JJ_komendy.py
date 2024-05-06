from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np
import argparse

o = object()

class Transformacje:
    def __init__(self, model: str = "WGS84"):
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
        if model == "WGS84":
            self.a = 6378137.0 # semimajor_axis
            self.b = 6356752.31424518 # semiminor_axis
        elif model == "GRS80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self.a = 6378245.0
            self.b = 6356863.01877  #tu można by dopisać po przecinku jeszcze liczby jak znajdziemy
        else:
            raise NotImplementedError(f"{model} model nie jest zaimplementowany")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2) # eccentricity  WGS84:0.0818191910428 
        self.ecc2 = (2 * self.flat - self.flat ** 2) # eccentricity**2
        
        
        """
        Poniższe funkcje są używane w celach pomocniczych dla obliczeń transformacji
        """
    def Npu(self, phi):
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
        return(N)

    def Sigma(self, phi):
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2)**2)/64) -  ((5*(self.ecc2)**3)/256)
        A2 = (3/8) * (self.ecc2 + ((self.ecc2)**2 / 4) + ((15*(self.ecc2)**3) / 128))
        A4 = (15/256) * ((self.ecc2)**2 + ((3*((self.ecc2)**3)) / 4))
        A6 = (35 * (self.ecc2)**3) / 3072
        sigma = self.a * ( A0 * phi - A2 * np.sin(2*phi) + A4 * np.sin(4*phi) - A6 * np.sin(6*phi))
        return(sigma)
        
        """
        Funkcje transforacji współrzędnych
        """
        
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
    
    
    
        # XYZ ---> NEU
        """
            Obliczenie macierzy Rneu
        """
    def Rneu(self, phi, lam):
        Rneu = np.array([[-np.sin(phi) * np.cos(lam), -np.sin(lam), np.cos(phi) * np.cos(lam)],
                         [-np.sin(phi) * np.sin(lam),  np.cos(lam), np.cos(phi) * np.sin(lam)],
                         [             np.cos(phi),            0,             np.sin(phi)]])
        return(Rneu)
    
    
        """
            Przeliczenie wsp XYZ na neu
        """
    def xyz2neu(self, X, Y, Z, X0, Y0, Z0):
        neu = []
        p = np.sqrt(X0**2 + Y0**2)
        phi = np.arctan(Z0 / (p * (1 - self.ecc2)))
        while True:
            N = self.Npu(phi)
            h = (p / np.cos(phi)) - N
            phi_poprzednia = phi
            phi = np.arctan((Z0 / p)/(1-((N * self.ecc2)/(N + h))))
            if abs(phi_poprzednia - phi) < (0.000001/206265):
                break 
        N = self.Npu(phi)
        h = p/np.cos(phi) - N
        lam = np.arctan(Y0 / X0)
        
        R_neu = self.Rneu(phi, lam)
        X_sr = [X - X0, Y - Y0, Z - Z0] 
        X_rneu = R_neu.T@X_sr
        neu.append(X_rneu.T)
            
        return(neu)



        # TRANSFORMACJA WSP BL ---> 1992
        """
            Algorytm przelicza współrzędne geodezyjne (BL) na współrzędne w układzie 1992 (XY)
        """
    def blto92(self, phi, lam):
        lam0 = (19 * np.pi)/180
        m = 0.9993
        wsp = []
        for phi,lam in zip(phi,lam):
            b2 = (self.a**2) * (1-self.ecc2)   #krotsza polowa
            e2p = (self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dlam = lam - lam0
            t = np.tan(phi)
            ni = np.sqrt(e2p * (np.cos(phi))**2)
            N = self.Npu(phi)
            sigma = self.Sigma(phi)
            
            xgk = sigma + ((dlam**2)/2) * N * np.sin(phi) * np.cos(phi) * ( 1+ ((dlam**2)/12)*(np.cos(phi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4))  + ((dlam**4)/360)*(np.cos(phi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2)))
            ygk = (dlam * N * np.cos(phi)) * (1+(((dlam)**2/6) * (np.cos(phi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(phi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)))
                        
            x92 = xgk * m - 5300000
            y92 = ygk * m + 500000
            wsp.append([x92, y92]) 
            
        return(wsp)
            
            
            
        # TRANSFORMACJA WSP BL ---> 2000
        """
            Następujący algorytm umożliwia przeliczenie współrzędnych geodezyjnych (BLH) na współrzędne w układzie 2000 (XY)
        """

    def blto00(self, phi, lam):
        m=0.999923
        print(phi, lam)
        wsp = []
        for phi,lam in zip(phi,lam):
            lam0=0 
            strefa = 0
            if lam >np.deg2rad(13.5) and lam < np.deg2rad(16.5):
                strefa = 5
                lam0 = np.deg2rad(15)
            elif lam >np.deg2rad(16.5) and lam < np.deg2rad(19.5):
                strefa = 6
                lam0 = np.deg2rad(18)
            elif lam >np.deg2rad(19.5) and lam < np.deg2rad(22.5):
                strefa =7
                lam0 = np.deg2rad(21)
            elif lam >np.deg2rad(22.5) and lam < np.deg2rad(25.5):
                strefa = 8
                lam0 = np.deg2rad(24)
            else:
                print("Punkt poza strefami odwzorowawczymi układu PL-2000")        
                     
            b2 = (self.a**2) * (1-self.ecc2)   #krotsza polos
            e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
            dlam = lam - lam0
            t = np.tan(phi)
            ni = np.sqrt(e2p * (np.cos(phi))**2)
            N = self.Npu(phi)
            sigma = self.Sigma(phi)
        
            xgk = sigma + ((dlam**2)/2) * N * np.sin(phi) * np.cos(phi) * ( 1+ ((dlam**2)/12)*(np.cos(phi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4))  + ((dlam**4)/360)*(np.cos(phi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2)))
            ygk = (dlam * N * np.cos(phi)) * (1+(((dlam)**2/6) * (np.cos(phi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(phi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)))
                     
            x00 = xgk * m
            y00 = ygk * m + strefa*1000000 + 500000
            wsp.append([x00, y00])
        return(wsp)  
        
    def pliczek(self, plik, funkcja):
        data = np.genfromtxt(plik,  delimiter = " ")
        if funkcja == "XYZ_BLH":
            X = data[:,0]
            Y = data[:,1]
            Z = data[:,2]
                # to zmienic e starej wersji, bedzie latwiej
            blh = self.hirvonen(X, Y, Z)
            np.savetxt(f"WYNIK_{funkcja}.txt", blh, delimiter=";")

        elif funkcja == "BLH_XYZ":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            h = data[:,2]
            XYZ = self.filh2XYZ(fi, lam, h)
            np.savetxt(f"WYNIK_{funkcja}.txt",XYZ, delimiter=";")
            
        elif funkcja == "XYZ_NEU":
            X0 = data[0,0]
            Y0 = data[0,1]
            Z0 = data[0,2]
            X = data[1,0]
            Y = data[1,1]
            Z = data[1,2]
            neu = self.xyz2neup(X, Y, Z, X0, Y0, Z0)
            np.savetxt(f"WYNIK_{funkcja}.txt", neu, delimiter=";")
            
        elif funkcja == "BL_PL1992":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            wsp92 = self.cale92(fi, lam)
            np.savetxt(f"WYNIK_{funkcja}.txt", wsp92, delimiter=";")
            
        elif funkcja == "BL_PL2000":
            fi = np.deg2rad(data[:,0])
            lam = np.deg2rad(data[:,1])
            wsp00 = self.cale00(fi, lam)
            np.savetxt(f"WYNIK_{funkcja}.txt", wsp00, delimiter=";")


# if __name__ == "__main__":
#     # utworzenie obiektu
#     geo = Transformacje(model = "WGS84")
#     print(sys.argv)
#     # dane XYZ geocentryczne
#     # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
#     # phi, lam, h = geo.xyz2plh(X, Y, Z)
#     # print(phi, lam, h)
#     # phi, lam, h = geo.xyz2plh2(X, Y, Z)
#     # print(phi, lam, h)
    

# with open ('wsp_inp.txt') as f:
# 	lines = f.readlines()
# 	lines = lines[4:]
# #	pprint(lines)

# 	coords_plh = []
# 	for line in lines:
# 		line = line.strip()
# 		x_str, y_str, z_str = line.split(',')
# 		x, y, z = (float(x_str), float(y_str), float(z_str))
# 		p, l, h = geo.xyz2plh(x, y, z)
# 		coords_plh.append([p, l, h])


# with open ('result_xyz2plh.txt', 'w') as f:
# 	for coords in coords_plh:
# 		coords_plh_line = ',' .join([str(coord) for coord in coords_plh])
# 		f.writelines(coords_plh_line + '\n')

# #dlaczego jest powielona jedna wartosc kilkukrotnie?


if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser(description="Podaj plik")
        parser.add_argument("-plik", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
        parser.add_argument("-elip", type = str, help = "Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'WGS84', 'GRS80', 'Elipsoida Krasowskiego' ")
        parser.add_argument("-funkcja", type = str, help = "Wybierz transformację jaką chcesz obliczyć: 'xyz2plh', 'plh2xyz', 'xyz2neu', 'blto92', 'blto00' ")
        args = parser.parse_args()
    except SyntaxError:
        print("Niestety nie ma takiego pliku. Spróbuj podać pełną scieżkę do pliku lub upewnij się że wpisujesz dobrą nazwę")
                   
    
    elip = {'WGS84':[6378137.000, 0.00669438002290], 'GRS80':[6378137.000, 0.00669438002290], 'Elipsoida Krasowskiego':[6378245.000, 0.00669342162296]}
    funkcja = {'XYZ_PLH' : 'xyz2plh', 'PLH_XYZ' : 'plh2xyz', 'XYZ_NEU' : 'xyz2neu', 'BL_PL1992' : 'blto92', 'BL_PL2000' : 'blto00'}
        
    try:
        geo = Transformacje(elip[args.elip.upper()])
        finito = geo.pliczek(args.plik, args.funkcja.upper())
        print("Zapisano")
    except KeyError:
        print("Podana funkcja/elipsoida nie istnieją, proszę upewnij się, że korzystasz z istniejących elipsoid")
    except AttributeError:
        print("Podana funkcja/elipsoida nie istnieje, proszę wprowadzić dostępne wartosci.")
    except FileNotFoundError:
        print("Nie znaleziono takiego pliku. Proszę spróbować wprowadzić inny plik.")
    except IndexError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")
    except ValueError:
        print("Nieodpowiedni format pliku. Proszę wprowadzić dane do pliku tak jak pokazano w przyładzie.")