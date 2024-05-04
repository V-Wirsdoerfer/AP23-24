import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp


### Daten generieren
Kennlinie_U, Kennlinie_pulses, Strom_Kennlinie = np.genfromtxt(
    "content/kennlinie.txt", unpack=True
)  # in V, 60, µA

# Zählraten

N_1 = (
    ufloat(154865, np.sqrt(154865)) / 120
)  # über 120s gemessen. deswegen aud s normieren
N_2 = (
    ufloat(145314, np.sqrt(145314)) / 120
)  # über 120s gemessen. deswegen aud s normieren
N_1_2 = (
    ufloat(258114, np.sqrt(258114)) / 120
)  # über 120s gemessen. deswegen aud s normieren

print("Die Messdaten für die Totzeit:")
print(N_1)
print(N_2)
print(N_1_2)


# Daten mit Fehler abspeichern
Kennlinie_pulses = unp.uarray(Kennlinie_pulses, np.sqrt(Kennlinie_pulses))


# in SI umrechnen
Strom_Kennlinie *= 1e-6  # in A

# Zählraten von pro 60s auf pro 1s umrechnen
Kennlinie_pulses = Kennlinie_pulses * 1 / 60


# Spannung gegen Strom plotten
fig, ax = plt.subplots()
ax.errorbar(
    Kennlinie_U,
    Strom_Kennlinie,
    yerr=0.05e-6,
    fmt="x",
    capsize=2,
    label="Messdaten",
)
ax.set(
    xlabel=r"Betriebsspannung $U$ / V",
    ylabel=r"Stromstärke $I$ / A",
)
ax.legend()
fig.savefig("build/Stromstaerke.pdf")


# Zählrate plotten
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
params_Kennlinie, cov_Kennlinie = np.polyfit(
    Kennlinie_U[4:19],
    unp.nominal_values(Kennlinie_pulses)[4:19],
    1,
    cov=True,
    w=[1 / i for i in unp.std_devs(Kennlinie_pulses)[4:19]],
)
ax.plot(
    Kennlinie_U[4:19],
    params_Kennlinie[0] * Kennlinie_U[4:19] + params_Kennlinie[1],
    label="Ausgleichsgerade",
)


ax.set(
    xlabel=r"Betriebsspannung $U$ / V",
    ylabel=r"Zählrate pro $S$",
    ylim=[15000 / 60, 19000 / 60],
    xlim=[300, 750],
)
ax.legend()
fig.savefig("build/Kennlinie.pdf")


### Steigung des Plateus berechnen mit Formel

print(
    Kennlinie_U[6],
    Kennlinie_U[9],
    Kennlinie_U[12],
    "\n",
    Kennlinie_pulses[6],
    Kennlinie_pulses[9],
    Kennlinie_pulses[12],
)
U60p = Kennlinie_pulses[12]
U60m = Kennlinie_pulses[6]
UA = Kennlinie_pulses[9]
s = (U60p - U60m) / (UA)
print("Die über die Formel berechnete Steigung des Plateaus ist: ", s)

Steigung_Kennlinie = ufloat(params_Kennlinie[0], np.sqrt(np.diag(cov_Kennlinie))[0])
print(
    "Die über die Ausgleichsgerade berechnete Steigung des Plateaus ist: ",
    Steigung_Kennlinie,
)

### Berechnung der Totzeit

t_tot = (N_1 + N_2 - N_1_2) / ((N_1_2) ** 2 - (N_1) ** 2 - (N_2) ** 2)
print("Das ist die Totzeit: ", t_tot, " s")
# print(N_1)
