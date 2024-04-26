from cProfile import label
from enum import auto
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp


### Konstanten
E0 = 4e6        #in eV
E0 *= 1.602e-19 #in J
x4 = 0.04       #in m
x5 = 0.05       #in m
p0 = 1.013      #in Bar | Normaldruck


### Daten generieren
counts_hist = np.genfromtxt("content/Verteilung.txt", unpack = True)
sum_pulses4, channel4, p4 = np.genfromtxt("content/name4.txt", unpack = True)   # in 1, 1, mBar, bei 4 cm
sum_pulses5, channel5, p5 = np.genfromtxt("content/name5.txt", unpack = True)   # in 1, 1, mBar, bei 5 cm


#in SI umrechnen
p4 *= 1e-3      #in Bar umrechnen
p5 *= 1e-3      #in Bar umrechnen



### Energie eines Channels ausrechnen
Betrag4 = channel4[0]
Betrag5 = channel5[0]


### Effektive Weglänge berechnen
eff_weglaenge4 = x4 * p4 / p0
eff_weglaenge5 = x5 * p5 / p0


### Aus Channel auf Energie Schließen
E4 = channel4 / Betrag4 * E0
E5 = channel5 / Betrag5 * E0


### Funktionen

def Steigung(params, cov):
    return ufloat(params[0], np.sqrt(abs(np.diag(cov)))[0])

def Achsenabschnitt(params, cov):
    return ufloat(params[1], np.sqrt(abs(np.diag(cov)))[1] )


### Druck gegen Energie auftragen
fig1,ax1 = plt.subplots(layout="constrained")
ax1.plot(eff_weglaenge4, E4,"x", label = "x = 4cm")
params4, cov4 = np.polyfit(eff_weglaenge4, E4, deg=1, cov=True)
x = np.linspace(0,0.025)
ax1.plot(x, unp.nominal_values(Steigung(params4, cov4) * x + Achsenabschnitt(params4, cov4)), label="Ausgleichsgerade")
ax1.legend()
ax1.set(
    xlabel=r"effektive Weglänge $x$",
    ylabel=r"Energie $E$"
)
fig1.savefig("build/4cm.pdf")

print("Die Ableitung dE/dx bei 4cm ist: ", Steigung(params4, cov4))


fig2,ax2 = plt.subplots(layout="constrained")
ax2.plot(eff_weglaenge5, E5,"x", label = "x = 5cm")
params5, cov5 = np.polyfit(eff_weglaenge5, E5, deg=1, cov=True)
ax2.plot(x, unp.nominal_values(Steigung(params5, cov5) * x + Achsenabschnitt(params5, cov5)), label="Ausgleichsgerade")
ax2.legend()
ax1.set(
    xlabel=r"effektive Weglänge $x$",
    ylabel=r"Energie $E$"
)
fig2.savefig("build/5cm.pdf")

print("Die Ableitung dE/dx bei 5cm ist: ", Steigung(params5, cov5))

### Zufallswerte erstellen
rng = np.random.default_rng(1)
mean = np.mean(counts_hist) #von wiki geklaut die Formel
poisson_Verteilung = rng.poisson(mean, size=100)
gauss_Verteilung = rng.normal(mean, 100, size=100)



#counts mit Histogramm plotten
fig3, axs = plt.subplots(2, 2, layout="constrained")
ax = np.ravel(axs)


#Histogramme in schleife erstellen
for i in np.arange(4):
    ax[i].hist(counts_hist, label="Messung", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].hist(poisson_Verteilung, label="Poisson", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].hist(gauss_Verteilung, label="Gauss", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].legend()
    ax[i].set_title(f"{10+2*i} - Bins")

fig3.savefig("build/Verteilung.pdf")

