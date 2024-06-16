import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

### Literaturwerte

e0   = 1.602 * 10**(-19)       # in Coulomb
m0   = 9.109 * 10**(-31)       # in kg
h    = 6.626 * 10**(-34)       # in Js
kB   = 1.381 * 10**(-34)       # in m²*kg/(s²*K)
eps0 = 8.854 * 10**(-12)       # in A*s*V*m

### Daten generieren

U_0, I_0_Anode  = np.genfromtxt("Kennlinie2.0.txt", unpack=True)
U_1, I_1_Anode  = np.genfromtxt("Kennlinie2.1.txt", unpack=True)
U_2, I_2_Anode  = np.genfromtxt("Kennlinie2.2.txt", unpack=True)
U_3, I_3_Anode  = np.genfromtxt("Kennlinie2.3.txt", unpack=True)
U_4, I_4_Anode  = np.genfromtxt("Kennlinie2.4.txt", unpack=True)

U_neg, I_Anl = np.genfromtxt("Anlaufstrom.txt", unpack=True)

### Daten in SI umrechen 

I_0_Anode   *= 1e-3
I_1_Anode   *= 1e-3
I_2_Anode   *= 1e-3
I_3_Anode   *= 1e-3
I_4_Anode   *= 1e-3
I_Anl       *= 1e-9

### Daten mit Fehlern 

U_0   = unp.uarray(U_0, 0.5)
U_1   = unp.uarray(U_1, 0.5)
U_2   = unp.uarray(U_2, 0.5)
U_3   = unp.uarray(U_3, 0.5)
U_4   = unp.uarray(U_4, 0.5)
U_neg = unp.uarray(U_neg, 0.5)

I_H0 = ufloat(2.0, 0.1)
I_H1 = ufloat(2.1, 0.1)
I_H2 = ufloat(2.2, 0.1)
I_H3 = ufloat(2.3, 0.1)
I_H4 = ufloat(2.4, 0.1)

### Kennlinien plotten

# Heizstrom von 2.0A

fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.errorbar(
    unp.nominal_values(U_0),
    I_0_Anode,
    xerr = unp.std_devs(U_0),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = "Messdaten mit Fehlerbalken",
)
ax1.set(
    xlabel = r"Spannung $U$ / V",
    ylabel = r"Strom $I$ / mA"
)
ax1.legend()
fig.savefig("build/Kennlinie_2.0.pdf")

# Heizstrom von 2.1A 

fig, ax2 = plt.subplots(1, 1, layout="constrained")
ax2.errorbar(
    unp.nominal_values(U_1),
    I_1_Anode,
    xerr = unp.std_devs(U_1),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = "Messdaten mit Fehlerbalken"
)
ax2.legend()
fig.savefig("build/Kennlinie_2.1.pdf")

# Heizstrom von 2.2A 

fig, ax3 = plt.subplots(1, 1, layout="constrained")
ax3.errorbar(
    unp.nominal_values(U_2),
    I_2_Anode,
    xerr = unp.std_devs(U_2),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = "Messdaten mit Fehlerbalken"
)
ax3.legend()
fig.savefig("build/Kennlinie_2.2.pdf")

# Heizstrom von 2.3A 

fig, ax4 = plt.subplots(1, 1, layout="constrained")
ax4.errorbar(
    unp.nominal_values(U_3),
    I_3_Anode,
    xerr = unp.std_devs(U_3),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = "Messdaten mit Fehlerbalken"
)
ax4.legend()
fig.savefig("build/Kennlinie_2.3.pdf")

# Heizstrom von 2.4A und Raumladungsgesetz

fig, ax5 = plt.subplots(1, 1, layout="constrained")
ax5.errorbar(
    unp.nominal_values(U_4),
    I_4_Anode,
    xerr = unp.std_devs(U_4),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = "Messdaten mit Fehlerbalken"
)

def I(U, a, b):
    return a * U ** b
U = np.linspace(0, unp.nominal_values(U_4[25]), 100)
params_Raumladung, cov_Raumladung = curve_fit(I, unp.nominal_values(U_4), I_4_Anode)
ax5.plot(U, I(U, *params_Raumladung))

ax5.legend()
fig.savefig("build/Kennlinie_2.4.pdf")
print(params_Raumladung[0], params_Raumladung[1])

### Anlaufstromgebiet

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(unp.nominal_values(U_neg), I_Anl, "x", label="Messdaten ohne Fehler")

def f(U, a, b):
    return a * np.exp(U * b)
U = np.linspace(min(unp.nominal_values(U_neg)), max(unp.nominal_values(U_neg)), 1000)
params_Anlauf, cov_Anlauf = curve_fit(f, unp.nominal_values(U_neg), I_Anl, p0=[4.5e-9, -1])
#ax.plot(U, f(U, *params_Anlauf))
print(params_Anlauf)
print(unp.nominal_values(U_neg).shape)
print(I_Anl.shape)
ax.legend()
fig.savefig("build/Anlaufstrom.pdf")


