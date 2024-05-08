import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

### Literaturwerte

c_Wasser = 1497        #m/s       
c_Acryl = 2730         #m/s 
c_Aluminium = 6320     #m/s

### Daten generieren

Hoehe_Alu, Laufzeit_Alu_2MHz = np.genfromtxt("content/Aluminium2MHz.txt", unpack=True)
Hoehe_Acryl_2MHz, Laufzeit_Acryl_2MHz, U_0_Acryl_2MHz, U_Acryl_2MHz = np.genfromtxt("content/Acryl2MHz.txt", unpack=True)
Hoehe_Acryl_1MHz, U_0_Acryl_1MHz, U_Acryl_1MHz = np.genfromtxt("content/Acryl1MHz.txt", unpack=True)
Fuellstand, Laufzeit_Fuell = np.genfromtxt("content/Kalibrierung.txt", unpack=True)
 
### Daten mit Fehlern 

#Hoehe_Alu = unp.uarray(Hoehe_Alu, 0.01)
#Hoehe_Acryl_1MHz = unp.uarray(Hoehe_Acryl_1MHz, 0.01)
#Hoehe_Acryl_2MHz = unp.uarray(Hoehe_Acryl_2MHz, 0.01)

### In SI umrechnen 

Hoehe_Alu *= 1e-3
Hoehe_Acryl_1MHz *= 1e-3
Hoehe_Acryl_2MHz *= 1e-3
Laufzeit_Acryl_2MHz *= 1e-6
Laufzeit_Alu_2MHz *= 1e-6

### Schallgeschwindigkeit bestimmen 

### Aluminium, 2MHz

fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(Laufzeit_Alu_2MHz, 2 * Hoehe_Alu, "rx", label="Messdaten")

#ax1.errorbar(
#    Laufzeit_Alu_2MHz,
#    unp.nominal_values(Hoehe_Alu),
#    yerr=unp.std_devs(Hoehe_Alu),
#    fmt="rx",
#    capsize=2,
#    label="Messdaten",
#)

params_WegZeitAlu, cov_WegZeitAlu = np.polyfit(Laufzeit_Alu_2MHz, 2 * Hoehe_Alu, deg=1, cov=True)
ax1.plot(Laufzeit_Alu_2MHz, params_WegZeitAlu[0] * Laufzeit_Alu_2MHz + params_WegZeitAlu[1], label="Ausgleichsgerade")

ax1.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $s$ / m"
)

ax1.legend()
fig.savefig("build/Schall_Alu.pdf")
print(params_WegZeitAlu[0]) 

### Acryl, 2MHz

fig, ax2 = plt.subplots(1, 1, layout="constrained")
ax2.plot(Laufzeit_Acryl_2MHz, 2 * Hoehe_Acryl_2MHz, "rx", label="Messdaten")
params_WegZeitAcryl, cov_WegZeitAcryl = np.polyfit(Laufzeit_Acryl_2MHz, 2 * Hoehe_Acryl_2MHz, deg=1, cov=True)
ax2.plot(Laufzeit_Acryl_2MHz, params_WegZeitAcryl[0] * Laufzeit_Acryl_2MHz + params_WegZeitAcryl[1], label="Ausgleichsgerade")

ax2.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $s$ / m"
)

ax2.legend()
fig.savefig("build/Schall_Acryl.pdf")
print(params_WegZeitAcryl[0])

### Kalibrierkurve erstellen 

fig, ax3 = plt.subplots(1, 1, layout="constrained")
ax3.plot(Fuellstand, Laufzeit_Fuell, "rx", label="Kalibrierdaten")
params_Kalibrier, cov_Kalibrier = np.polyfit(Fuellstand, Laufzeit_Fuell, deg=2, cov=True)
ax3.plot(Fuellstand, params_Kalibrier[0] * Fuellstand**2  + params_Kalibrier[1] * Fuellstand + params_Kalibrier[2], label="Kalibrierkurve")
ax3.legend()
fig.savefig("build/Kalibrierkurve.pdf")

print(params_Kalibrier)
