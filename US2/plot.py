import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat


###Grundgrößen und Konstanten
c_Acryl = 2730      #in m/s
c_Wasser = 1497     #in m/s
c_Linse = 2500      #in m/s
c_Glaskörper = 1410 #in m/s
H_Acrylblock = ufloat(0.079,0.000025)


### Daten generieren
H_top, H_bot = np.genfromtxt("./content/acrylblock.txt", unpack=True)         #in mm und mm
t_top, t_bot = np.genfromtxt("./content/Laufzeit_Acryl.txt", unpack=True)     #in µs und µs
Herz_T, Herz_Amp = np.genfromtxt("./content/Herz.txt", unpack=True)           #in  s und µs


### Daten in SI umrechnen
H_top *= 1e-3       #von mm in m
H_bot *= 1e-3       #von mm in m

t_top *= 1e-6       #von µs in s
t_bot *= 1e-6       #von µs in s

Herz_Amp *= 1e-6    #von µs in s



### Durchmesser mithilfe der Gemessenen Höhen berechnen
Durchmesser_Messschieber = H_Acrylblock - H_top - H_bot
print("Die Durchmesser über den Messschieber bestimmt: ", Durchmesser_Messschieber)

### Durchmesser mithilfe der Laufzeiten berechnen
Durchmesser_Laufzeiten = H_Acrylblock - t_top * c_Acryl - t_bot * c_Acryl
print("Die Durchmesser über die Laufzeiten bestimmt: ", Durchmesser_Laufzeiten)


### Abweichungen bestimmen
Durchmesser_Differenz = Durchmesser_Laufzeiten - Durchmesser_Messschieber
print("Die Differenz zwischen den Methoden der Durchmesser ist: ", Durchmesser_Differenz)



