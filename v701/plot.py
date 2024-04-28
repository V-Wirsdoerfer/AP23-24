from cProfile import label
from enum import auto
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp


### Konstanten
E0 = 4e6        #in eV
#E0 *= 1.602e-19 #in J
x4 = 0.04       #in m
x5 = 0.05       #in m
p0 = 1.013      #in Bar | Normaldruck
p0 *= 1e5       #von Bar in Pascal


### Daten generieren
counts_hist = np.genfromtxt("content/Verteilung.txt", unpack = True)
sum_pulses4, channel4, p4 = np.genfromtxt("content/name4.txt", unpack = True)   # in 1, 1, mBar, bei 4 cm
sum_pulses5, channel5, p5 = np.genfromtxt("content/name5.txt", unpack = True)   # in 1, 1, mBar, bei 5 cm


#in SI umrechnen
p4 *= 1e-3      #in Bar umrechnen
p5 *= 1e-3      #in Bar umrechnen

p4 *= 1e5       #in Pascal umrechnen
p5 *= 1e5       #in Pascal umrechnen



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
ax1.plot(eff_weglaenge4, E4,"x", label = "Messdaten")
params4, cov4 = np.polyfit(eff_weglaenge4, E4, deg=1, cov=True)
x = np.linspace(0,0.025)
ax1.plot(x, unp.nominal_values(Steigung(params4, cov4) * x + Achsenabschnitt(params4, cov4)), label="Ausgleichsgerade")
ax1.legend()
ax1.set(
    xlabel=r"effektive Weglänge $x$ / m",
    ylabel=r"Energie $E$ / eV"
)
fig1.savefig("build/4cm.pdf")

print("Die Ableitung dE/dx bei 4cm ist: ", Steigung(params4, cov4))


fig2,ax2 = plt.subplots(layout="constrained")
ax2.plot(eff_weglaenge5, E5,"x", label = "Messdaten")
params5, cov5 = np.polyfit(eff_weglaenge5, E5, deg=1, cov=True)
ax2.plot(x, unp.nominal_values(Steigung(params5, cov5) * x + Achsenabschnitt(params5, cov5)), label="Ausgleichsgerade")
ax2.legend()
ax2.set(
    xlabel=r"effektive Weglänge $x$ / m",
    ylabel=r"Energie $E$ / eV"
)
fig2.savefig("build/5cm.pdf")

print("Die Ableitung dE/dx bei 5cm ist: ", Steigung(params5, cov5))

### Zufallswerte erstellen
rng = np.random.default_rng(1)
mean = np.mean(counts_hist) #von wiki geklaut die Formel
poisson_Verteilung = rng.poisson(mean, size=100)
gauss_Verteilung = rng.normal(mean, 100, size=100)

print("Mittelwert/lambda = ", mean)
print("Varianz: ", np.std(counts_hist)**2)

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



### Counts gegen freie Weglänge auftragen
## 4cn
fig4, ax4 = plt.subplots(layout="constrained")
ax4.plot(eff_weglaenge4, sum_pulses4, "x", label="Messdaten")
ax4.hlines((max(sum_pulses4) - min(sum_pulses4)) / 2 + min(sum_pulses4) , 0, 0.025, label="Hälfte der gemessenen Pulses")
ax4.set(
    xlabel = r"effektive Weglänge $x$ / m",
    ylabel = "Nummer an counts", 
)
ax4.legend(loc="lower left")
fig4.savefig("build/sumpulses4.pdf")

##5cm
fig5, ax5 = plt.subplots(layout="constrained")
ax5.plot(eff_weglaenge5, sum_pulses5, "x", label="Messdaten")
ax5.hlines((max(sum_pulses5) - min(sum_pulses5)) / 2 + min(sum_pulses5), 0, 0.025, label="Hälfte der gemessenen Pulses")
ax5.set(
    xlabel = r"effektive Weglänge $x$ / m",
    ylabel = "Nummer an counts", 
)
ax5.legend(loc="lower left")
fig5.savefig("build/sumpulses5.pdf")


