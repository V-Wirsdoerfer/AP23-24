import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

### Daten generieren/Konstanten

z_Wellenlaenge, z_EVK, z_BLF = np.genfromtxt("Zaehlung.txt", unpack=True, dtype=int)
p0 = ufloat(1010e2, 0.01e-2)        # in pascal
T = ufloat(23.7, 0.05) + 273.15     # in Kelvin
lambda_lit = 635e-9                 # in m
D = 50e-3                           # in m
#delta_p = 500                       # in mmHg
delta_p = 66661.195                 # in Pascal

### Wellenlänge berechnen

x = ufloat(5e-3, 0.01e-3)   # in m

z_Wellenlaenge_korr = z_Wellenlaenge[z_Wellenlaenge!=z_Wellenlaenge[2]]
z_Wellenlaenge_mean = ufloat(np.mean(z_Wellenlaenge_korr), np.std(z_Wellenlaenge_korr))
Wellenlaenge = (2 * x) / (5.017 * z_Wellenlaenge_mean)
print("Durchschnittliche Zählrate: ", z_Wellenlaenge_mean)
print("Die Wellenlänge beträgt: ", Wellenlaenge, " m")

#Abweichung berechnen
delta_Wellenlaenge = abs(635e-9 - Wellenlaenge)/635e-9
print("Die berechnete Wellenlänge weicht um ", delta_Wellenlaenge*100, "% ab.")



### Berechnung des Brechungsindex
z_Brechungsindex = np.concatenate(( (z_EVK[0:5]), (z_BLF[0:5]) ))
z_Brechungsindex_mean = ufloat(np.mean(z_Brechungsindex), np.std(z_Brechungsindex))
delta_n = z_Brechungsindex_mean * lambda_lit / (2* D)
n = 1 + delta_n * T * p0 / (273.15 * delta_p)

print("Der Brechungsindex von Luft beträgt: ", n)

#Abweichung berechnen
Abweichung_n = abs(1.00029-n)/1.00029
print("Die prozentuale Abweichung der Brechungsindizes zum Literaturwert beträgt: ", Abweichung_n*100, "%.")