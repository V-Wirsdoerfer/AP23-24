import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.optimize import curve_fit

### Literaturwerte

e0   = 1.602 * 10**(-19)       # in Coulomb
m0   = 9.109 * 10**(-31)       # in kg
h    = 6.626 * 10**(-34)       # in Js
kB   = 1.381 * 10**(-23)       # in m²*kg/(s²*K)
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

U_H0 = ufloat(4, 0.5)
U_H1 = ufloat(4, 0.5)
U_H2 = ufloat(4.5, 0.5)
U_H3 = ufloat(5, 0.5)
U_H4 = ufloat(5, 0.5)


### Kennlinien plotten

# Heizstrom von 2.0 - 2.3 A

fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.errorbar(
    unp.nominal_values(U_0),
    I_0_Anode,
    xerr = unp.std_devs(U_0),
    yerr = None,
    fmt = "r.", 
    capsize= 0, 
    label = r"$I_H = 2.0\,$A"
)

ax1.errorbar(
    unp.nominal_values(U_1),
    I_1_Anode,
    xerr = unp.std_devs(U_1),
    yerr = None,
    fmt = "b.",
    capsize= 0,
    label = r"$I_H = 2.1\,$A"
)
ax1.errorbar(
    unp.nominal_values(U_2),
    I_2_Anode,
    xerr = unp.std_devs(U_2),
    yerr = None,
    fmt = "y.",
    capsize= 0,
    label = r"$I_H = 2.2\,$A"
)

ax1.errorbar(
    unp.nominal_values(U_3),
    I_3_Anode,
    xerr = unp.std_devs(U_3),
    yerr = None,
    fmt = "g.",
    capsize= 0,
    label = r"$I_H = 2.3\,$A"
)

ax1.set(
    xlabel = r"Spannung $U$ / V",
    ylabel = r"Strom $I$ / mA"
)
ax1.legend()
fig.savefig("build/Kennlinie_2.0.pdf")

# Heizstrom von 2.4 A

fig, ax5 = plt.subplots(1, 1, layout="constrained")
ax5.errorbar(
    unp.nominal_values(U_4),
    I_4_Anode,
    xerr = unp.std_devs(U_4),
    yerr = None,
    fmt = "r.",
    capsize= 0,
    label = r"$I_H = 2.4\,$A"
)

ax5.legend()
fig.savefig("build/Kennlinie_2.4.pdf")

### Raumladungsgesetz über den Logarithmus 

log_I_4_Anode = np.log(I_4_Anode[1:25])
log_U_4 = np.log(unp.nominal_values(U_4)[1:25])

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(log_U_4, log_I_4_Anode, "rx", label="logarithmierte Messdaten")

params_Raumladung, cov_Raumladung = np.polyfit(
    log_U_4, 
    log_I_4_Anode, 
    deg=1, 
    cov=True,
 )
ax.plot(log_U_4, params_Raumladung[0] * log_U_4 + params_Raumladung[1])

ax.legend()
fig.savefig("build/logRaum.pdf")

m_log_Raum = ufloat(params_Raumladung[0], np.sqrt(np.diag(cov_Raumladung))[0])
b_log_Raum = ufloat(params_Raumladung[1], np.sqrt(np.diag(cov_Raumladung))[1])

print("Parameter logarithmiertes Raumladungsgesetz:\n")
print("Steigung der Ausgleichsgerade:\n", m_log_Raum)
print("Achsenabschnitt der Ausgleichgerade:\n", b_log_Raum)


### Anlaufstromgebiet

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(unp.nominal_values(U_neg), I_Anl, "x", label="Messdaten ohne Fehler")

def f(U, a, b):
    return a * np.exp(U * b)
U = np.linspace(min(unp.nominal_values(U_neg)), max(unp.nominal_values(U_neg)), 1000)
params_Anlauf, cov_Anlauf = curve_fit(f, unp.nominal_values(U_neg), I_Anl, p0=[5e-9, -4])
ax.plot(U, f(U, *params_Anlauf))

ax.legend()
fig.savefig("build/Anlaufstrom.pdf")

#print("Achsenabschnitt exponentieller fit:\n", params_Anlauf[0])
#print("Exponentialkoeffizient des fits:\n",params_Anlauf[1])


### Logarithmiertes Anlaufgebiet

log_U_neg = np.log(unp.nominal_values(U_neg)[1:-1])
log_I_Anl = np.log(I_Anl[1:-1])

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(unp.nominal_values(U_neg)[1:-1], log_I_Anl, "rx", label="logarithmierte Messdaten")

params_AnlaufPolyfit, cov_AnlaufPolyfit = np.polyfit(
    unp.nominal_values(U_neg)[1:-1], 
    log_I_Anl, 
    deg=1, 
    cov=True,
 )
ax.plot(unp.nominal_values(U_neg)[1:-1], params_AnlaufPolyfit[0] * unp.nominal_values(U_neg)[1:-1] + params_AnlaufPolyfit[1])

ax.legend()
fig.savefig("build/logAnlauf.pdf")

m_log_Anl = ufloat(params_Anlauf[0], np.sqrt(np.diag(cov_Anlauf))[0])
b_log_Anl = ufloat(params_Anlauf[1], np.sqrt(np.diag(cov_Anlauf))[1])

print("Parameter logarithmiertes Anlaufgebiet:\n")
print("Steigung der Gerade:\n", m_log_Anl)
print("Achsenabschnitt der Gerade:\n", b_log_Anl)

#T = (-e0) / (kB * params_AnlaufPolyfit[0])
#print(T)


### Berechnung der Temperatur
N_WL = ufloat(0.95,0.05) #zwischen 0.9 und 1    #in W
f_Diode = 0.35      #in cm^2
eta = 0.28          # Emissionsgrad
sigma = 5.7e-12     #in W/(cm^2 K^4)

def Temperatur(I, U):
    T = ( (I*U - N_WL)/(f_Diode * eta * sigma) )**(1/4)
    return T

U_H = [U_H0, U_H1, U_H2, U_H3, U_H4]
I_H = [I_H0, I_H1, I_H2, I_H3, I_H4]

for i in range(len(U_H)):
    print("Die Temperatur für eine Heizspannung von ", U_H[i], "V ist ", Temperatur(I_H[i], U_H[i]))
#for (U, I) in zip((), ()):
 #   print("U:\n", U)
