import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 

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

### Fehler einbeziehen 

H_top = unp.uarray(H_top, 0.025e-3)
H_bot = unp.uarray(H_bot, 0.025e-3)

### Anpassungsschicht berechnen

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(t_top, 2 * unp.nominal_values(H_top), "rx", label="Messdaten")
params_top, cov_top = np.polyfit(t_top, 2 * unp.nominal_values(H_top), deg=2, cov=True)
ax.plot(t_top, params_top[0] * t_top + params_top[1], label="Ausgleichsgerade")
ax.legend()
plt.show()

### Durchmesser mithilfe der Gemessenen Höhen berechnen
Durchmesser_Messschieber = H_Acrylblock - H_top - H_bot
print("Die Durchmesser über den Messschieber bestimmt: ", Durchmesser_Messschieber)

### Durchmesser mithilfe der Laufzeiten berechnen
Durchmesser_Laufzeiten = H_Acrylblock - t_top * c_Acryl - t_bot * c_Acryl
print("Die Durchmesser über die Laufzeiten bestimmt: ", Durchmesser_Laufzeiten)


### Abweichungen bestimmen
Durchmesser_Differenz = Durchmesser_Laufzeiten - Durchmesser_Messschieber
print("Die Differenz zwischen den Methoden der Durchmesser ist: ", Durchmesser_Differenz)



