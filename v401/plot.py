import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

### Daten generieren/Konstanten

z_Wellenlaenge, z_EVK, z_BLF = np.genfromtxt("Zaehlung.txt", unpack=True)
p0 = ufloat(1010e2, 0.01e-2)         # in pascal
T0 = ufloat(23.7, 0.05) + 273.15   # in Kelvin

### Wellenlänge berechnen

x = ufloat(5e-3, 0.01e-3)   # in m

z_Wellenlaenge_korr = z_Wellenlaenge[z_Wellenlaenge!=z_Wellenlaenge[2]]
z_Wellenlange_mean = ufloat(np.mean(z_Wellenlaenge_korr), np.std(z_Wellenlaenge_korr))
Wellenlange = (2 * x) / (5.017 * z_Wellenlange_mean)
print("Durchschnittliche Zählrate: ", z_Wellenlange_mean)
print("Die Wellenlänge beträgt: ", Wellenlange, " m")

#Abweichung berechnen
delta_Wellenlaenge = abs(635e-9 - Wellenlange)/635e-9
print("Die berechnete Wellenlänge weicht um ", delta_Wellenlaenge*100, "% ab.")

### Berechnung des Brechungsindex


