from math import sin, cos, sqrt, atan, atan2, degrees, radians
import sys
import numpy as np
import argparse

o = object()

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
        
    def xyz2plh(self, X, Y, Z):
        wyniki = []
        for X, Y, Z in zip(X, Y, Z):
            p = np.sqrt(X**2+Y**2)
            phi = np.arctan(Z/(p*(1-self.ecc2)))
            while True:
                N = self.Npu(phi)
                h = (p/np.cos(phi)) - N
                phi_prev = phi
                phi = np.arctan((Z / p)/(1 - ((N * self.ecc2) / (N + h))))
                if abs(phi_prev - phi)<(0.000001/206265):
                    break
            N = self.Npu(phi)
            h = p / np.cos(phi) - N
            lam = np.arctan(Y / X)
            wyniki.append([np.rad2deg(phi), np.rad2deg(lam), h])
        return wyniki
    
    
    def plh2xyz(self, phi, lam, h):
        phi = radians(phi)
        lam = radians(lam)
        Rn = self.a / sqrt(1 - self.ecc2 * sin(phi)**2)
        q = Rn * self.ecc2 * sin(phi)
        x = (Rn + h) * cos(phi) * cos(lam)
        y = (Rn + h) * cos(phi) * sin(lam)
        z = (Rn + h) * sin(phi) * q
        return x, y, z
    
        
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
        wyniki = []
        for phi, lam in zip (phi, lam):   
            b2 = self.a**2 * (1 - self.ecc2)
            ep2 = (self.a**2 - b2) / b2
            dlam = lam - lam0
            t = np.tan(phi)
            ni2 = ep2 * (np.cos(phi)**2)
            N = self.Npu(phi)
         
            A0 = 1- (self.ecc2 / 4) - (3 * self.e2**2 / 64) - (5 * self.e2**3 / 256)
            A2 = (3/8) * (self.ecc2 + (self.ecc2**2 / 4) + (15 * self.ecc2**3 / 128))
            A4 = (15/256) * (self.ecc2**2 + ((3 * self.e2**3) / 4))
            A6 = (35 * self.ecc2**3) / 3072
            sigma = self.a * (A0 * phi - A2 * np.sin(2*phi) + A4 * np.sin(4 * phi) - A6 * np.sin( 6 * phi))
        
            xgk =  sigma + (((dlam**2 / 2) * N * np.sin(phi) * np.cos(phi)) * (1 + ((dlam**2 / 12) * (np.cos(phi)**2) * (5 - t**2 + 9 * ni2 + 4 * ni2**2)) + ((dlam**4 / 360) * (np.cos(phi)**4) * (61 - 58 * t**2 + t**4 + 270 * ni2 - 330 * ni2 * t**2))))
            ygk =  (dellam* N * np.cos(fi))  *   ( 1 +  ((dellama**2/6)   *   (np.cos(fi)**2)   *  (1 - t**2 + ni2))     +     (((dellama**4/120)*(np.cos(fi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
        
            x92 = xgk * m - 5300000
            y92 = ygk*m + 500000
            wyniki.append([x92,y92])

        return  wyniki

    def BLto2000(self,fi,lama,m=0.999923):
        wyniki = []
        for fi, lama in zip (fi,lama):
            lama0 = 0
            strefa = 0
            if lama >np.deg2rad(13.5) and lama < np.deg2rad(16.5):
                strefa = 5
                lama0 = np.deg2rad(15)
            elif lama >np.deg2rad(16.5) and lama < np.deg2rad(19.5):
                strefa = 6
                lama0 = np.deg2rad(18)
            elif lama >np.deg2rad(19.5) and lama < np.deg2rad(22.5):
                strefa =7
                lama0 = np.deg2rad(21)
            elif lama >np.deg2rad(22.5) and lama < np.deg2rad(25.5):
                strefa = 8
                lama0 = np.deg2rad(24)
            
            b2 = self.a**2*(1-self.e2)    
            ep2 = (self.a**2-b2)/b2
            dellama = lama - lama0
            t = np.tan(fi)
            ni2 = ep2*(np.cos(fi)**2)
            N = self.Npu(fi)
             
            A0 = 1- (self.e2/4)-(3*self.e2**2/64)-(5*self.e2**3/256)
            A2 = (3/8)*(self.e2+(self.e2**2/4)+(15*self.e2**3/128))
            A4 = (15/256)*(self.e2**2+((3*self.e2**3)/4))
            A6 = (35*self.e2**3)/3072
            
            sigma = self.a *(A0*fi-A2*np.sin(2*fi)+A4*np.sin(4*fi)-A6*np.sin(6*fi))
            
            xgk =  sigma    +    ( ((dellama**2/2)*N*np.sin(fi)*np.cos(fi))    *    (1   +   ((dellama**2/12)*(np.cos(fi)**2)*(5 - t**2 + 9*ni2 + 4*ni2**2))      +         ((dellama**4/360)*(np.cos(fi)**4)*(61 - 58*t**2 + t**4 + 270*ni2 - 330*ni2*t**2))))
            ygk =  (dellama* N * np.cos(fi))  *   ( 1 +  ((dellama**2/6)   *   (np.cos(fi)**2)   *  (1 - t**2 + ni2))     +     (((dellama**4/120)*(np.cos(fi)**4)) * (5 - (18*t**2) + t**4 + (14 * ni2) - (58*ni2*t**2))))
            
            x2000 = xgk * m 
            y2000 = ygk*m + (strefa *1000000) +500000
            wyniki.append([x2000,y2000])
            
        return  wyniki  
        

if __name__ == "__main__":
    # utworzenie obiektu
    geo = Transformacje(model = "WGS84")
    # dane XYZ geocentryczne
    # X = 3664940.500; Y = 1409153.590; Z = 5009571.170
    # phi, lam, h = geo.xyz2plh(X, Y, Z)
    # print(phi, lam, h)
    # phi, lam, h = geo.xyz2plh2(X, Y, Z)
    # print(phi, lam, h)
    
    input_file_path = sys.argv[-1]
    if '--header_lines' in sys.argv:
        header_lines = int(sys.argv[2])
    if '--xyz2plh' in sys.argv and '--plh2xyz' in sys.argv:
        print('Możesz podać tylko jedną flagę')
        
    elif '--xyz2plh' in sys.argv:
        with open (input_file_path, 'r') as f:
         	lines = f.readlines()
         	lines = lines[header_lines]


         	coords_plh = []
         	for line in lines:
                 line = line.strip()
                 x_str, y_str, z_str = line.split(',')
                 x, y, z = (float(x_str), float(y_str), float(z_str))
                 p, l, h = geo.xyz2plh(x, y, z)
                 coords_plh.append([p, l, h])

        with open ('result_xyz2plh.txt', 'w') as f:
            f.write('phi[deg], lam[deg], h[m] \n')
            for coords in coords_plh:
                coords_plh_line = ','.join([str(coord) for coord in coords])
                f.writelines(coords_plh_line + '\n')
                
    elif '--plh2xyz' in sys.argv:
        with open (input_file_path, 'r') as f:
         	lines = f.readlines()
         	lines = lines[header_lines]

         	coords_xyz = []
         	for line in lines:
                 line = line.strip()
                 phi_str, lam_str, h_str = line.split(',')
                 phi, lam, h = (float(phi_str), float(lam_str), float(h_str))
                 x, y, z = geo.xyz2plh(phi, lam, h)
                 coords_xyz.append([x, y, z])

        with open ('result_plh2xyz.txt', 'w') as f:
            f.write('x[m], y[m], z[m] \n')
            for coords in coords_xyz:
                coords_xyz_line = ','.join([f'{coord:11.3f}' for coord in coords])
                f.write(coords_xyz_line + '\n')
                
                
    elif '--xyz2neu' in sys.argv:
        with open (input_file_path, 'r') as f:
         	lines = f.readlines()
         	lines = lines[header_lines]


         	coords_neu = []
         	for line in lines:
                 line = line.strip()
                 x, y, z = line.split(',')
                 x, y, z = (float(x), float(y), float(z))
                 x_0, y_0, z_0 = [float(coord) for coord in sys.argv[-4:-1]]
                 n, e, u = geo.xyz2neu(x, y, z, x_0, y_0, z_0)
                 coords_neu.append([n, e, u])

        with open ('result_xyz2neu.txt', 'w') as f:
            f.write('n[m], e[m], u[m] \n')
            for coords in coords_neu:
                coords_neu_line = ','.join([f'{coord:11.3f}' for coord in coords])
                f.writelines(coords_neu_line + '\n')



