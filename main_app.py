# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 12:04:30 2022

@author: Julia
"""

from __future__ import unicode_literals
import sys
from PyQt5.QtWidgets import QDialog, QApplication
from Aplikacja_geodezyjna import *
from math import sin, cos, tan, sqrt, pi


class MyApp(QDialog):
    def __init__(self):
        super().__init__()
        self.ui = Ui_Aplikacja_geodezyjna()
        self.ui.setupUi(self)
        self.ui.GRS80.clicked.connect(self.get_radiobutton)
        self.ui.WGS84.clicked.connect(self.get_radiobutton)
        self.ui.Strefa_5.clicked.connect(self.strefa)
        self.ui.Strefa_6.clicked.connect(self.strefa)
        self.ui.Strefa_7.clicked.connect(self.strefa)
        self.ui.Strefa_8.clicked.connect(self.strefa)
        self.ui.xy92.clicked.connect(self.fl2xy92)
        self.ui.xy2000.clicked.connect(self.fl2xy00)
        self.ui.xy_gk.clicked.connect(self.fl2xy_gk00)
        self.setWindowIcon(QtGui.QIcon('ikona.jpg'))
        self.show()
        

# ----------------------------------------------------------------------------------------------------------------------  

    def get_radiobutton(self):
        """
        Funkcja przyjmująca parametry elipsoidy odnieseinia GRS84 lub WGS80 
        w zależnosci od wyboru użytkownika w oknie dialogowym.

        """
        radioButton = self.sender()
        if radioButton.isChecked():
            self.ui.model1.setText("Model " + radioButton.text())
            self.ui.model2.setText("Model " + radioButton.text())
            self.ui.model3.setText("Model " + radioButton.text())
        if self.ui.GRS80.isChecked():
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif self.ui.WGS84.isChecked():
            self.a = 6378137.0
            self.b = 6356752.31424518
            
        self.flat = (self.a - self.b) / self.a
        self.ecc = sqrt(2 * self.flat - self.flat ** 2)
        self.ecc2 = (2 * self.flat - self.flat ** 2)
        self.b2 = self.b ** 2

        
# ---------------------------------------------------------------------------------------------------------------------  

    def strefa(self):
        """
        Funkcja przyjmująca wartosć południka srodkowego l0 na podstawie wyboru 
        użytkownika w oknie dialogowym.

        """
        radioButton = self.sender()
        if radioButton.isChecked():
            self.ui.Wyswietlenie_wyniku_xy00.setText("Strefa " + radioButton.text())
        if self.ui.Strefa_5.isChecked():
            self.l0 = 15*pi/180
        if self.ui.Strefa_6.isChecked():
            self.l0 = 18*pi/180
        if self.ui.Strefa_7.isChecked():
            self.l0 = 21*pi/180
        if self.ui.Strefa_8.isChecked():
            self.l0 = 24*pi/180
        
        
# ----------------------------------------------------------------------------------------------------------------------      
 
    def func_n(self, f: float) -> float:
        """
        Funkcja obliczająca promień krzywizny w I wertykale
 
        Parameters
        ----------                
        INPUT:
            phi : [float] - szerokość geodezyjna (decimal degrees)

        Returns
        -------     
        OUTPUT:
            N : [float] - promien krzywizny (meters)
        """
        N = (self.a)/sqrt(1-self.ecc2*(sin(f)**2))
        
        return(N)
    
# --------------------------------------------------------------------------------------------------------------    

    def sigma (self, f: float) -> float:        
        """
        Algorytm obliczający długosć łuku południka na podstawie podanej wartosci szerokosci geograficznej

        Parameters
        ----------
        INPUT:
            f : [float] - szerokość geodezyjna (decimal degrees)

        Returns
        -------            
        OUTPUT:
            si : [float] - długosć łuku południka (meters)
            
        """               
        A0 = 1 - (self.ecc2/4)-((3*(self.ecc2**2))/64)-((5*(self.ecc2**3))/256)
        A2 = (3/8)*(self.ecc2+((self.ecc2**2)/4)+((15*(self.ecc2**3))/128))
        A4 = (15/256)*((self.ecc2**2)+((3*(self.ecc2**3))/4))
        A6 = (35*(self.ecc2**3))/3072
        si = self.a*(A0*f - A2*sin(2*f) + A4*sin(4*f) - A6*sin(6*f))

        return(si)

# -------------------------------------------------------------------------------------------------------------------- 
    
    def fl2xy00(self):
        """
        Algorytm przeliczający współrzędne geodezyjne (phi, lambda) na współrzędne 
        w układzie PL-2000.
        
        Parameters
        ----------            
        INPUT:
            f : [float] - szerokość geodezyjna (decimal degrees)
            l : [float] - długość geodezyjna (decimal degreees)

        Returns
        -------    
        OUTPUT:
            x00 : [float] - współrzędna X w układzie PL-2000 (meters)
            y00 : [float] - współrzędna Y w układzie PL-2000 (meters)

        """
        if len(self.ui.Wartosc_phi.text())!=0:
            f = float(self.ui.Wartosc_phi.text())
        else:
            f = 0
        if len(self.ui.Wartosc_lambdy.text())!=0:
            l = float(self.ui.Wartosc_lambdy.text())
        else:
            l = 0

        f = f*pi/180
        l = l*pi/180
        ep2 = (self.a**2-self.b2)/self.b2
        t = tan(f)
        n2 = ep2*(cos(f)**2)
        N = self.func_n(f)
        si = self.sigma(f)
        dL = l - self.l0 
        xgk = si + ((dL**2)/2)*N*sin(f)*cos(f)*(1 + (dL**2/12)*cos(f)**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dL*N*cos(f)*(1 + (dL**2/6)*cos(f)**2*(1 - t**2 + n2) + (dL**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))            
        m2000 = 0.999923
        x00 = round(xgk * m2000, 3)
        y00 = round(ygk * m2000 + (self.l0*180/pi/3)* 1000000 + 500000, 3)
        self.ui.Wyswietlenie_wyniku_xy00.setText("x00: " + str(x00) + " [m]" + " ;  y00: " + str(y00) + " [m]")
 
# -----------------------------------------------------------------------------------------------------------------------    

    def fl2xy92(self):
        """
        Algorytm przeliczający współrzędne geodezyjne (phi, lambda) na współrzędne 
        w układzie PL-1992.

        Parameters
        ----------      
        INPUT:
            f : [float] - szerokość geodezyjna (decimal degreees)
            l : [float] - długość geodezyjna (decimal degreees)
    
        Returns
        -------    
        OUTPUT:
            x92 : [float] - współrzędna X w układzie PL-1992 (meters)
            y92 : [float] - współrzędna Y w układzie PL-1992 (meters)
    
        """
        if len(self.ui.Wartosc_phi.text())!=0:
            f = float(self.ui.Wartosc_phi.text())
        else:
            f = 0
        if len(self.ui.Wartosc_lambdy.text())!=0:
            l = float(self.ui.Wartosc_lambdy.text())
        else:
            l = 0
        
        f = f*pi/180
        l = l*pi/180
        ep2 = (self.a**2-self.b2)/self.b2
        t = tan(f)
        n2 = ep2*(cos(f)**2)
        N = self.func_n(f)
        si = self.sigma(f)    
        dL = l - 19*pi/180
        xgk = si + ((dL**2)/2)*N*sin(f)*cos(f)*(1 + (dL**2/12)*cos(f)**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2))
        ygk = dL*N*cos(f)*(1 + (dL**2/6)*cos(f)**2*(1 - t**2 + n2) + (dL**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2))
        x92 = round(xgk * 0.9993-5300000, 3)
        y92 = round(ygk *0.9993+500000, 3)
        self.ui.Wyswietlenie_wyniku_xy92.setText("x92: " + str(x92) + " [m]" + " ;  y92: " + str(y92) + " [m]" )
        
 # ------------------------------------------------------------------------------------------------------------------
    
    def get_l0(self, l: float) -> int:
         
        """
        Algorytm obliczający wartosć południka srodkowego L0 na podstawie l
        
        INPUT:
            l : [float] - długość geodezyjna (decimal degrees)
            
        OUTPUT:
            L0 : [float] - południk srodkowy w danym układzie (radians)
        """        
        if 13.5 < l <= 16.5:
            L0 = 15 * pi / 180
        if 16.5 < l <= 19.5:
            L0 = 18 * pi / 180
        if 19.5 < l <= 22.5:
            L0 = 21 * pi / 180
        if 22.5 < l <= 25.5:
            L0 = 24 * pi / 180
            
        return(L0)    
   
# --------------------------------------------------------------------------------------------------------------  

    def fl2xy_gk00(self):
        """
        Algorytm przeliczający współrzędne godezyjne (phi, lambda) na współrzędne 
        w odwzorowaniu Gaussa-Krugera (xgk, ygk)

        Parameters
        ----------
        INPUT:
            phi : [float] - szerokość geodezyjna (decimal degrees)
            lam : [float] - długość geodezyjna (decimal degrees)

        Returns
        -------
        OUTPUT:
            xgk :[float] - współrzędna x w odwzorowaniu Gaussa-Krugera (meters)
            ygk :[float] - współrzędna y w odwzorowaniu Gaussa-Krugera (meters)

        """
        if len(self.ui.Wartosc_phi.text())!=0:
            f = float(self.ui.Wartosc_phi.text())
        else:
            f = 0
        if len(self.ui.Wartosc_lambdy.text())!=0:
            l = float(self.ui.Wartosc_lambdy.text())
        else:
            l = 0

        f = f*pi/180
        l = l*pi/180
        ep2 = (self.a**2-self.b2)/self.b2
        t = tan(f)
        n2 = ep2*(cos(f)**2)
        N = self.func_n(f)
        si = self.sigma(f)
        dL = l - self.get_l0(l*180/pi) 
        xgk = round(si + ((dL**2)/2)*N*sin(f)*cos(f)*(1 + (dL**2/12)*cos(f)**2*(5 - t**2 + 9*n2 + 4*n2**2) + (dL**4/360)*cos(f)**4*(61 - 58*t**2 + t**4 + 270*n2 - 330*n2*t**2)),3)
        ygk = round(dL*N*cos(f)*(1 + (dL**2/6)*cos(f)**2*(1 - t**2 + n2) + (dL**4/120)*cos(f)**4*(5 - 18*t**2 + t**4 + 14*n2 - 58*n2*t**2)),3)
       
        self.ui.Wyswietlenie_wyniku_xy_gk.setText("xgk: " + str(xgk) + " [m]" + " ;  ygk: " + str(ygk) + " [m]")
                       
# ----------------------------------------------------------------------------------------------------------------------  

if __name__=="__main__":
	app = QApplication(sys.argv)
	w = MyApp()
	w.show()
	sys.exit(app.exec_())
