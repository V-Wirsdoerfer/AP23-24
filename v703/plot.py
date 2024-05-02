import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp



### Daten generieren
Kennlinie_U, Kennlinie_pulses, Strom_Kennlinie = np.genfromtxt("content/kennlinie.txt", unpack=True) # in V, 1, µA



#in SI umrechnen
Strom_Kennlinie *= 1e-6 #in A



#Spannung gegen Strom plotten
fig, ax = plt.subplots()
ax.errorbar(
    Kennlinie_U,
    Strom_Kennlinie,
    yerr = 0.05e-6,
    fmt = "x",
    capsize=2,
    label="Messdaten",
)

ax.set(
    xlabel=r"Betriebsspannung $U$ / V",
    ylabel =r"Stromstärke $I$ / A",
)

ax.legend()
fig.savefig("build/Stromstaerke.pdf")



fig, ax = plt.subplots()
ax.errorbar(
    Kennlinie_U,
    Kennlinie_pulses,
    yerr=np.sqrt(Kennlinie_pulses),
    fmt="rx",
    capsize=2,
    label="Messdaten",
)

### Ausgleichsgerade
params_Kennlinie, cov_Kennlinie = np.polyfit(Kennlinie_U[4:19], Kennlinie_pulses[4:19], 1, cov=True, w= [1 / i for i in np.sqrt(Kennlinie_pulses)[4:19]])
ax.plot(Kennlinie_U[4:19], params_Kennlinie[0] * Kennlinie_U[4:19] + params_Kennlinie[1], label="Ausgleichsgerade")


ax.set(
    xlabel=r"Betriebsspannung $U$ / V",
    ylabel="Zählrate",
    ylim=[15000, 19000],
    xlim=[300, 750]
)
ax.legend()
fig.savefig("build/Kennlinie.pdf")



### Steigung des Plateus berechnen mit Formel

print(
    Kennlinie_U[6],
    Kennlinie_U[9],
    Kennlinie_U[12],    
)
s = (Kennlinie_pulses[12]- Kennlinie_pulses[6]) / Kennlinie_pulses[9]

print("Die über die Formel berechnete Steigung des Plateaus ist: ", s)
Steigung_Kennlinie = ufloat(params_Kennlinie[0], np.sqrt(np.diag(cov_Kennlinie))[0])
print("Die über die Ausgleichsgerade berechnete Steigung des Plateaus ist: ", Steigung_Kennlinie)