import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

Kennlinie_U, val_Kennlinie_pulses, Strom_Kennlinie = np.genfromtxt("content/kennlinie.txt", unpack=True) # in V, 1, µA

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


### Fehler behaften
Kennlinie_pulses = unp.uarray(val_Kennlinie_pulses, np.sqrt(val_Kennlinie_pulses))

fig, ax = plt.subplots()
ax.errorbar(
    Kennlinie_U,
    unp.nominal_values(Kennlinie_pulses),
    yerr=unp.std_devs(Kennlinie_pulses),
    fmt="rx",
    capsize=2,
    label="Messdaten",
)

### Ausgleichsgerade

params_Kennlinie, cov_Kennlinie = np.polyfit(Kennlinie_U[4:19], val_Kennlinie_pulses[4:19], 1, cov=True, w= [1 / i for i in np.sqrt(val_Kennlinie_pulses)[4:19]])
ax.plot(Kennlinie_U[4:19], params_Kennlinie[0] * Kennlinie_U[4:19] + params_Kennlinie[1], label="Ausgleichsgerade")


ax.set(
    xlabel=r"Betriebsspannung $U$ / V",
    ylabel="Zählrate",
    ylim=[15000, 19000],
    xlim=[300, 750]
)
ax.legend()
fig.savefig("build/Kennlinie.pdf")
