from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np
import argparse


class Transformacje:
    def __init__(self, model: str = "WGS84"):
        if model == "WGS84":
            self.a = 6378137.0 
            self.b = 6356752.31424518 
        elif model == "GRS80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self.a = 6378245.0
            self.b = 6356863.01877  
        else:
            raise NotImplementedError(f"{model} model nie jest zaimplementowany")
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)  
        self.ecc2 = (2 * self.flat - self.flat ** 2) 
        
        
    def Npu(self, phi):
        phi = radians(phi)
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
        return(N)
        
        """
        Funkcje transforacji współrzędnych
        """
        
    def xyz2plh(self, X, Y, Z, output = 'dec_degree'):
        r   = sqrt(X**2 + Y**2)           # promień
        lat_prev = atan(Z / (r * (1 - self.ecc2)))    # pierwsze przybliilizenie
        lat = 0
        while abs(lat_prev - lat) > 0.000001/206265:    
            lat_prev = lat
            N = self.a / sqrt(1 - self.ecc2 * sin(lat_prev)**2)
            h = r / cos(lat_prev) - N
            lat = atan((Z/r) * (((1 - self.ecc2 * N/(N + h))**(-1))))
        lon = atan(Y/X)
        N = self.a / sqrt(1 - self.ecc2 * (sin(lat))**2);
        h = r / cos(lat) - N       
        if output == "dec_degree":
            return degrees(lat), degrees(lon), h 
        elif output == "dms":
            lat = self.deg2dms(degrees(lat))
            lon = self.deg2dms(degrees(lon))
            return f"{lat[0]:02d}:{lat[1]:02d}:{lat[2]:.2f}", f"{lon[0]:02d}:{lon[1]:02d}:{lon[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - output format not defined")
    
    
    def plh2xyz(self, phi, lam, h):
        phi = radians(phi)
        lam = radians(lam)
        X = (self.Npu(phi) + h) * cos(phi) * cos(lam)
        Y = (self.Npu(phi) + h) * cos(phi) * sin(lam)
        Z = (self.Npu(phi) * (1 - self.ecc2) + h) * sin(phi)
        return X, Y, Z
    
        
    def xyz2neu(self, x, y, z, x_0, y_0, z_0):
        phi, lam, _ = [radians(coord) for coord in self.xyz2plh(x, y, z)]
        
        R = np.array([[-sin(lam),  -sin(phi)*cos(lam), cos(phi)*cos(lam)],
                      [ cos(lam),  -sin(phi)*sin(lam), cos(phi)*sin(lam)],
                      [    0,          cos(phi),          sin(phi)     ]])
        
        xyz_t = np.array([[x - x_0],
                          [y - y_0],
                          [z - z_0]])
        
        [[E], [N], [U]] = R.T @ xyz_t
        
        return N, E, U



        # TRANSFORMACJA WSP BL ---> 1992
        """
            Algorytm przelicza współrzędne geodezyjne (BL) na współrzędne w układzie 1992 (XY)
        """
    def BLto92(self, phi, lam, m=0.9993):  
        lam0 = np.deg2rad(19)
        phi = np.radians(phi)
        lam = np.radians(lam)
        wyniki = []
        # for phi, lam in zip (phi, lam):   
        b2 = self.a**2 * (1 - self.ecc2)
        ep2 = (self.a**2 - b2) / b2
        dlam = lam - lam0
        t = np.tan(phi)
        ni2 = ep2 * (np.cos(phi)**2)
        N = self.Npu(phi)
         
        A0 = 1- (self.ecc2 / 4) - (3 * self.ecc2**2 / 64) - (5 * self.ecc2**3 / 256)
        A2 = (3/8) * (self.ecc2 + (self.ecc2**2 / 4) + (15 * self.ecc2**3 / 128))
        A4 = (15/256) * (self.ecc2**2 + ((3 * self.ecc2**3) / 4))
        A6 = (35 * self.ecc2**3) / 3072
        sigma = self.a * (A0 * phi - A2 * np.sin(2*phi) + A4 * np.sin(4 * phi) - A6 * np.sin( 6 * phi))
    
        xgk =  sigma + (((dlam**2 / 2) * N * np.sin(phi) * np.cos(phi)) * (1 + ((dlam**2 / 12) * (np.cos(phi)**2) * (5 - t**2 + 9 * ni2 + 4 * ni2**2)) + ((dlam**4 / 360) * (np.cos(phi)**4) * (61 - 58 * t**2 + t**4 + 270 * ni2 - 330 * ni2 * t**2))))
        ygk =  (dlam* N * np.cos(phi))  *   ( 1 +  ((dlam**2/6)   *   (np.cos(phi)**2)   *  (1 - t**2 + ni2))     +     (((dlam**4/120)*(np.cos(phi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
        
        x92 = xgk * m - 5300000
        y92 = ygk*m + 500000
        wyniki.append([x92,y92])

        return  wyniki

    def BLto2000(self,phi,lam,m=0.999923):
        wyniki = []
        phi = np.radians(phi)
        lam = np.radians(lam)
        lama0 = 0
        strefa = 0
        if lam >np.deg2rad(13.5) and lam < np.deg2rad(16.5):
            strefa = 5
            lama0 = np.deg2rad(15)
        elif lam >np.deg2rad(16.5) and lam < np.deg2rad(19.5):
            strefa = 6
            lama0 = np.deg2rad(18)
        elif lam >np.deg2rad(19.5) and lam < np.deg2rad(22.5):
            strefa =7
            lama0 = np.deg2rad(21)
        elif lam >np.deg2rad(22.5) and lam < np.deg2rad(25.5):
            strefa = 8
            lama0 = np.deg2rad(24)
            
        b2 = self.a**2 * (1 - self.ecc2)
        ep2 = (self.a**2 - b2) / b2
        dlam = lam - lama0
        t = np.tan(phi)
        ni2 = ep2 * (np.cos(phi)**2)
        N = self.Npu(phi)
        
        A0 = 1- (self.ecc2 / 4) - (3 * self.ecc2**2 / 64) - (5 * self.ecc2**3 / 256)
        A2 = (3/8) * (self.ecc2 + (self.ecc2**2 / 4) + (15 * self.ecc2**3 / 128))
        A4 = (15/256) * (self.ecc2**2 + ((3 * self.ecc2**3) / 4))
        A6 = (35 * self.ecc2**3) / 3072
        sigma = self.a * (A0 * phi - A2 * np.sin(2*phi) + A4 * np.sin(4 * phi) - A6 * np.sin( 6 * phi))
    
        xgk =  sigma + (((dlam**2 / 2) * N * np.sin(phi) * np.cos(phi)) * (1 + ((dlam**2 / 12) * (np.cos(phi)**2) * (5 - t**2 + 9 * ni2 + 4 * ni2**2)) + ((dlam**4 / 360) * (np.cos(phi)**4) * (61 - 58 * t**2 + t**4 + 270 * ni2 - 330 * ni2 * t**2))))
        ygk =  (dlam* N * np.cos(phi))  *   ( 1 +  ((dlam**2/6)   *   (np.cos(phi)**2)   *  (1 - t**2 + ni2))     +     (((dlam**4/120)*(np.cos(phi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
    
        x2000 = xgk * m 
        y2000 = ygk*m + (strefa *1000000) +500000
        wyniki.append([x2000,y2000])
            
        return  wyniki  
        

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "WGS84")
    g80 = Transformacje(model = "GRS80")
    kra = Transformacje(model = "Elipsoida Krasowskiego")
    try:
        if '--header_lines' not in sys.argv:
            print('BŁĄD: Nie podano flagi --header_lines')            
        else:
            input_file_path = sys.argv[-1]
            if '--header_lines' in sys.argv:
                header_lines = int(sys.argv[2])
            else:
                print("Podaj flagę --header_lines oraz od której linijki w pliku zaczynają się dane")
            
            dozwolone_flagi = ["--xyz2plh_geo", "--plh2xyz_geo", "--xyz2neu_geo", "--blto92_geo", "--blto2000_geo", "--xyz2plh_g80", "--plh2xyz_g80", 
                               "--xyz2neu_g80", "--blto92_g80", "--blto2000_g80","--xyz2plh_kra", "--plh2xyz_kra", "--xyz2neu_kra", "--blto92_kra", "--blto2000_kra"]
            if len([arg for arg in sys.argv if arg in dozwolone_flagi]) > 1:
                print("BŁĄD: Wpisano więcej niż jedną flagę.")

            elif '--xyz2plh_geo' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
                   
               	coords_plhgeo = []
               	for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = (float(x_str), float(y_str), float(z_str))
                    p, l, h = geo.xyz2plh(x, y, z)
                    coords_plhgeo.append([p, l, h])

                with open('result_xyz2plh_wgs84.txt', 'w') as f:
                    f.write('phi[deg], lam[deg], h[m] \n')
                    for coords in coords_plhgeo:
                        coords_plh_linegeo = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_plh_linegeo + '\n')

            elif '--plh2xyz_geo' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]

               	coords_xyzgeo = []
               	for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = (float(phi_str), float(
                        lam_str), float(h_str))
                    x, y, z = geo.xyz2plh(phi, lam, h)
                    coords_xyzgeo.append([x, y, z])

                with open('result_plh2xyz_wgs84.txt', 'w') as f:
                    f.write('x[m], y[m], z[m] \n')
                    for coords in coords_xyzgeo:
                        coords_xyz_linegeo = ','.join(
                            [str(coord)for coord in coords])
                        f.write(coords_xyz_linegeo + '\n')

            elif '--xyz2neu_geo' in sys.argv:
                coords_neugeo = []
                with open(input_file_path, 'r') as f:
              	  lines = f.readlines()
               	lines = lines[header_lines:]

               	for line in lines:
                    line = line.strip()
                    x, y, z = line.split(',')
                    x, y, z = (float(x), float(y), float(z))
                    x_0, y_0, z_0 = [float(coord)
                                     for coord in sys.argv[-4:-1]]
                    n, e, u = geo.xyz2neu(x, y, z, x_0, y_0, z_0)
                    coords_neugeo.append([n, e, u])

                with open('result_xyz2neu_wgs84.txt', 'w') as f:
                    f.write('n[m], e[m], u[m] \n')
                    for coords in coords_neugeo:
                        coords_neu_linegeo = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_neu_linegeo + '\n')

            elif '--blto92_geo' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl92geo = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = geo.BLto92(phi, lam)
                    x92, y92 = result[0]
                    coords_bl92geo.append([x92, y92])

                with open('result_blto92_wgs84.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl92geo:
                        coords_bl92_linegeo = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl92_linegeo + '\n')

            elif '--blto2000_geo' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl2000geo = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = geo.BLto2000(phi, lam)
                    x2000, y2000 = result[0]
                    coords_bl2000geo.append([x2000, y2000])

                with open('result_blto2000_wgs84.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl2000geo:
                        coords_bl2000_linegeo = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl2000_linegeo + '\n')

        #elipsoida grs80
            elif '--xyz2plh_g80' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_plhg80 = []
               	for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = (float(x_str), float(y_str), float(z_str))
                    p, l, h = g80.xyz2plh(x, y, z)
                    coords_plhg80.append([p, l, h])

                with open('result_xyz2plh_grs80.txt', 'w') as f:
                    f.write('phi[deg], lam[deg], h[m] \n')
                    for coords in coords_plhg80:
                        coords_plh_lineg80 = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_plh_lineg80 + '\n')

            elif '--plh2xyz_g80' in sys.argv:
                with open(input_file_path, 'r') as f:
                    lines = f.readlines()
               	lines = lines[header_lines:]

               	coords_xyzg80 = []
               	for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = (float(phi_str), float(
                        lam_str), float(h_str))
                    x, y, z = g80.xyz2plh(phi, lam, h)
                    coords_xyzg80.append([x, y, z])

                with open('result_plh2xyz_grs80.txt', 'w') as f:
                    f.write('x[m], y[m], z[m] \n')
                    for coords in coords_xyzg80:
                        coords_xyz_lineg80 = ','.join(
                            [str(coord)for coord in coords])
                        f.write(coords_xyz_lineg80 + '\n')

            elif '--xyz2neu_g80' in sys.argv:
                coords_neug80 = []
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]

               	for line in lines:
                    line = line.strip()
                    x, y, z = line.split(',')
                    x, y, z = (float(x), float(y), float(z))
                    x_0, y_0, z_0 = [float(coord)
                                     for coord in sys.argv[-4:-1]]
                    n, e, u = g80.xyz2neu(x, y, z, x_0, y_0, z_0)
                    coords_neug80.append([n, e, u])

                with open('result_xyz2neu_grs80.txt', 'w') as f:
                    f.write('n[m], e[m], u[m] \n')
                    for coords in coords_neug80:
                        coords_neu_lineg80 = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_neu_lineg80 + '\n')

            elif '--blto92_g80' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl92g80 = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = g80.BLto92(phi, lam)
                    x92, y92 = result[0]
                    coords_bl92g80.append([x92, y92])

                with open('result_blto92_gr80.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl92g80:
                        coords_bl92_lineg80 = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl92_lineg80 + '\n')

            elif '--blto2000_g80' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl2000g80 = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = g80.BLto2000(phi, lam)
                    x2000, y2000 = result[0]
                    coords_bl2000g80.append([x2000, y2000])

                with open('result_blto2000_grs80.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl2000g80:
                        coords_bl2000_lineg80 = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl2000_lineg80 + '\n')

        #elipsoida krasowskiego
            elif '--xyz2plh_kra' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_plhkra = []
               	for line in lines:
                    line = line.strip()
                    x_str, y_str, z_str = line.split(',')
                    x, y, z = (float(x_str), float(y_str), float(z_str))
                    p, l, h = kra.xyz2plh(x, y, z)
                    coords_plhkra.append([p, l, h])

                with open('result_xyz2plh_kras.txt', 'w') as f:
                    f.write('phi[deg], lam[deg], h[m] \n')
                    for coords in coords_plhkra:
                        coords_plh_linekra = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_plh_linekra + '\n')

            elif '--plh2xyz_kra' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]

               	coords_xyzkra = []
               	for line in lines:
                    line = line.strip()
                    phi_str, lam_str, h_str = line.split(',')
                    phi, lam, h = (float(phi_str), float(
                        lam_str), float(h_str))
                    x, y, z = kra.xyz2plh(phi, lam, h)
                    coords_xyzkra.append([x, y, z])

                with open('result_plh2xyz_kras.txt', 'w') as f:
                    f.write('x[m], y[m], z[m] \n')
                    for coords in coords_xyzkra:
                        coords_xyz_linekra = ','.join(
                            [str(coord)for coord in coords])
                        f.write(coords_xyz_linekra + '\n')

            elif '--xyz2neu_kra' in sys.argv:
                coords_neukra = []
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]

               	for line in lines:
                    line = line.strip()
                    x, y, z = line.split(',')
                    x, y, z = (float(x), float(y), float(z))
                    x_0, y_0, z_0 = [float(coord)
                                     for coord in sys.argv[-4:-1]]
                    n, e, u = kra.xyz2neu(x, y, z, x_0, y_0, z_0)
                    coords_neukra.append([n, e, u])

                with open('result_xyz2neu_kras.txt', 'w') as f:
                    f.write('n[m], e[m], u[m] \n')
                    for coords in coords_neukra:
                        coords_neu_linekra = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_neu_linekra + '\n')

            elif '--blto92_kra' in sys.argv:
                with open(input_file_path, 'r') as f:
                 	lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl92kra = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = kra.BLto92(phi, lam)
                    x92, y92 = result[0]
                    coords_bl92kra.append([x92, y92])

                with open('result_blto92_kras.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl92kra:
                        coords_bl92_linekra = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl92_linekra + '\n')

            elif '--blto2000_kra' in sys.argv:
                with open(input_file_path, 'r') as f:
               	  lines = f.readlines()
               	lines = lines[header_lines:]
               	coords_bl2000kra = []

               	for line in lines:
                    line = line.strip()
                    phi, lam, h = line.split(',')
                    phi, lam, h = (float(phi), float(lam), float(h))
                    result = kra.BLto2000(phi, lam)
                    x2000, y2000 = result[0]
                    coords_bl2000kra.append([x2000, y2000])

                with open('result_blto2000_kras.txt', 'w') as f:
                    f.write('x[m], y[m] \n')
                    for coords in coords_bl2000kra:
                        coords_bl2000_linekra = ','.join(
                            [str(coord) for coord in coords])
                        f.writelines(coords_bl2000_linekra + '\n')

    except IndexError:
        print('BŁĄD: nie podano odpowiednich wartosci')
    except SyntaxError:
        print('BŁĄD: nie podano odpowiednich wartosci')
    except FileNotFoundError:
        print('BŁĄD: Niepoprawnie podana nazwa pliku z danymi początkowymi')
    except ValueError:
        print('Podaj współrzędne geocentryczne srodka układu (neu)')
    finally:
        print('Brawo, program wykonał zadanie poprawnie a plik został utworzony pod nazwą "result_<nazwa transformacji".txt')
 


