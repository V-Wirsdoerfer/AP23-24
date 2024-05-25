import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

### Daten generieren/Konstanten

z_Wellenlaenge, z_EVK, z_BLF = np.genfromtxt("Zaehlung.txt", unpack=True)
p0 = 1010e2,         # in pascal
T0 = 23.7 + 273.15   # in Kelvin

### Wellenlänge berechnen

x = ufloat(5e-3, 0.01e-3)   # in m

z_Wellenlaenge_korr = z_Wellenlaenge[z_Wellenlaenge!=z_Wellenlaenge[2]]
z_Wellenlange_mean = np.mean(z_Wellenlaenge_korr)
Wellenlange = (2 * x) / (5.017 * z_Wellenlange_mean)
print("Die Wellenlänge beträgt: ", Wellenlange, " m")

### Berechnung des Brechungsindex


