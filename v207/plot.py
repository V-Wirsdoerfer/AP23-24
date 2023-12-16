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


# Funktionen für statiische Viskosität, Apparaturkonstante und dynamische Viskosität
def get_eta_kl(t):
    return K_kl * (rho_kl - rho_Fl) * t


def get_K_gr(t):
    return eta_kl.nominal_value / ((rho_gr - rho_Fl) * t)


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
    xlabel="Fallzeiten",
    ylabel="Viskosität",
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
    xlabel="Fallzeiten",
    ylabel="Viskosität",
)
ax2.legend()
fig2.savefig("build/gross_20C.pdf")


# plt.show()
