from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
import scipy
from uncertainties import ufloat
from scipy.stats import sem

# from scipy.optimize import curve_fitpyplot as plt

# Daten generieren
tg_ou1, tg_ou2, tg_uo1, tg_uo2 = np.genfromtxt("./content/gross_20C.txt", unpack="True")
tk_ou, tk_uo = np.genfromtxt("./content/klein_20C.txt", unpack="True")
T, td_ou1, td_ou2, td_uo1, td_uo2 = np.genfromtxt(
    "./content/dynamisch.txt", unpack="True"
)

# Daten: Statisch klein plotten
K_kl = 0.07640
rho_kl = (4.9528) / (4 / 3 * np.pi * (1.576 / 2) ** 3)
rho_gr = 4.4531 / (4 / 3 * np.pi * (0.5 * 1.599) ** 3)
rho_Fl = 0.99821

print("Dichte große Kugel = ", rho_gr)
print("Dichte kleine Kugel = ", rho_kl)


# Funktionen für statiische Viskosität, Apparaturkonstante und dynamische Viskosität
def get_eta_kl(t):
    return K_kl * (rho_kl - rho_Fl) * t


def get_eta_gr(t):
    return K_gr * (rho_gr - rho_Fl) * t


def get_K_gr(t):
    return eta_kl.nominal_value / (
        (rho_gr - rho_Fl) * t
    )  # evtl noch t -> 2t, da nur die Hälfte der Zeit


######################################
# Die Auswertung erfolgt der Reihenfolge nach. Erst werden die Daten der kleinen Kugel gesichtet,
# dann wird aus der kleinen Kugel die Viskosität über einen Mittwelwert aller Viskositäten berechnet.
# Anschließend werden die berechneten Viskositäten gesichtet. Dann wird aus diesen die Apparaturkonstante
# Für diegroße Kugel berechnet und gemittelt.
######################################

# Kontrollplot 20 grad kleine Kugel
fig1, ax1 = plt.subplots(
    layout="constrained", label="Kontrollplot kleine Kugel 20 grad"
)
ax1.plot(tk_ou, get_eta_kl(tk_ou), "rx", label="oben nach unten")
ax1.plot(tk_uo, get_eta_kl(tk_uo), "gx", label="unten nach oben")

ax1.set(
    xlabel=r"$t / \unit{s}$",
    ylabel=r"$\eta / \unit{\gram \per \centi \meter \per \second}$",
)
ax1.legend()
fig1.savefig("build/klein_20C.pdf")

# Viskositäten berechnen und dann über alle Mitteln
ar_eta_kl = np.concatenate((np.array(get_eta_kl(tk_ou)), np.array(get_eta_kl(tk_uo))))
eta_kl = ufloat(np.mean(ar_eta_kl), scipy.stats.sem(ar_eta_kl))
print("Das ist die Viskosität über die kleine Kugel berechnet: ", eta_kl)


# Apparaturkonstante aus den einzelnen Viskositäten berechnen
ar_K_gr = np.concatenate(
    (get_K_gr(tg_ou1), get_K_gr(tg_ou2), get_K_gr(tg_uo1), get_K_gr(tg_uo2))
)
K_gr = ufloat(np.mean(ar_K_gr), scipy.stats.sem(ar_K_gr))
print("Das ist die Apparaturkonstante der großen Kugel: ", K_gr)

# Kontrollplot 20 grad große Kugel
fig2, ax2 = plt.subplots(label="Kontrollplot große Kugel 20 grad")
ax2.plot(tg_ou1, get_K_gr(tg_ou1), "rx", label="1 oben nach unten")
ax2.plot(tg_ou2, get_K_gr(tg_ou2), "gx", label="2 oben nach unten")
ax2.plot(tg_uo1, get_K_gr(tg_uo1), "bx", label="1 unten nach oben")
ax2.plot(tg_uo2, get_K_gr(tg_uo2), "kx", label="2 unten nach oben")

ax2.set(
    xlabel=r"$t / s$",
    ylabel=r"$K / \unit{\milli\pascal\centi\cubic\meter\per\gram}$",
)
ax2.legend()
fig2.savefig("build/gross_20C.pdf")


# dynamische Viskosität


# Daten linearisieren
def lin_eta(t):
    eta = K_gr.nominal_value * (rho_gr - rho_Fl) * t
    return np.log(eta)


def lin_T(T):
    return 1 / T


fig3, ax3 = plt.subplots(label="dynamische Viskosität")
ax3.plot(lin_T(T), lin_eta(td_ou1), "rx", label="1 oben nach unten")
ax3.plot(lin_T(T), lin_eta(td_ou2), "gx", label="2 oben nach unten")
ax3.plot(lin_T(T), lin_eta(td_uo1), "bx", label="1 unten nach oben")
ax3.plot(lin_T(T), lin_eta(td_uo2), "kx", label="2 unten nach oben")


# Ausgleichsgerade berechnen, indem erst gemittelt wird.
average_t = np.zeros(np.size(T))
for i in range(np.size(T)):
    average_t[i] = np.mean((td_ou1[i], td_ou2[i], td_uo1[i], td_uo2[i]))
    # print(average_t)
# print(average_t)
# print(lin_eta(average_t))
params, cov = np.polyfit(lin_T(T), lin_eta(average_t), deg=1, cov=True)
x = np.linspace(0.018, 0.044)
ax3.plot(x, params[0] * x + params[1], label="Ausgleichsgerade")
ax3.set(
    xlabel=r"$\frac{1}{T} / \unit{\per\second}$",
    ylabel=r"$\ln{\eta}$",
    xlim=(0.018, 0.044),
)
ax3.legend()
fig3.savefig("build/dynamisch.pdf")
cov = abs(np.diag(cov))
err = np.sqrt(
    cov
)  # Cov darf warum auch immer nicht direkt in der sqrt berechnet werden. war pain das herauszufinden (╥_╥)
ln_A = ufloat(params[1], err[1])
A = np.e**ln_A
B = ufloat(params[0], err[0])
print("Koeffizient A = ", A, "\nKoeffizient B = ", B)


# Reynoldzahl berechnen
ar_tg = np.concatenate((tg_ou1, tg_ou2, tg_uo1, tg_uo2))
tg = ufloat(np.mean(ar_tg), scipy.stats.sem(ar_tg))
Re_g = rho_Fl * 5 / tg * 1.576 / eta_kl
print(
    "\nDaten zum Rechnen, große Kugel:",
    "\nDichte: ",
    rho_Fl,
    "\nFallstrecke: ",
    5,
    "\nFallzeit: ",
    tg,
    "\nFallgeschwindigkeit: ",
    5/tg,
    "\nLänge L: ",
    1.576,
    "\nViskosität des Wassers:",
    eta_kl,
)


ar_tk = np.concatenate((tk_ou, tk_uo))
tk = ufloat(np.mean(ar_tk), scipy.stats.sem(ar_tk))
Re_k = (rho_Fl * (10 / tk) * 1.559) / eta_kl
print(
    "\nDaten zum Rechnen, kleine Kugel:",
    "\nDichte: ",
    rho_Fl,
    "\nFallstrecke: ",
    10,
    "\nFallzeit: ",
    tk,
    "\nFallgeschwindigkeit:",
    10/tk,
    "\nLänge L: ",
    1.559,
    "\nViskosität des Wassers:",
    eta_kl,
)
print("\nReynold klein: ", Re_k, "\nReynold groß: ", Re_g)


tg = ufloat(np.mean(ar_tg), scipy.stats.sem(ar_tg))
Re_g = rho_Fl * 5 / (average_t) * 1.576 / get_eta_gr(average_t)
print(
    "\nDaten zum Rechnen, dynamisch große Kugel:",
    "\nDichte: ",
    rho_Fl,
    "\nFallstrecke: ",
    5,
    "\nFallzeit: ",
    average_t,
    "\nFallgeschwindigkeit:",
    5/average_t,
    "\nLänge L: ",
    1.576,
    "\nViskosität des Wassers:",
    get_eta_gr(average_t),
)
print("\nReynoldzahl dynamisch: ", Re_g)


# plt.show()
