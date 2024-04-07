from cProfile import label
from distutils.command import sdist
from turtle import distance
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from uncertainties import ufloat
import uncertainties.unumpy as unp


# Funktionen zur Auswertung
# Fehlerrechnung
def mittelwert(x):
    return sum(x) / np.size(x)


def Fehler_mittelwert(x):
    x_M = mittelwert(x)
    sum = 0
    for i in range(np.size(x)):
        sum += (x[i] - x_M) ** 2
    x_F = np.sqrt(sum / (np.size(x) * (np.size(x) - 1)))
    return x_F


def kappa(rho, c, A_nah, A_fern, dt_a):
    dx = ufloat(0.03, 0.005)
    A_nah = ufloat(mittelwert(A_nah), Fehler_mittelwert(A_nah))
    A_fern = ufloat(mittelwert(A_fern), Fehler_mittelwert(A_fern))
    dt = ufloat(mittelwert(dt_a), Fehler_mittelwert(dt_a))
    return (rho * c * (dx**2)) / (2 * dt * unp.log(A_nah / A_fern))


def Waermestrom(kappa, A, T):
    dx = ufloat(0.03, 0.005)
    return -kappa * A * (T / dx)


# Ausgabefunktion
def Wert(s, x):
    print(s, ": ")
    print("Mittelwert von ", s, ": ", mittelwert(x))
    print("Fehler des Mittwelwertes von ", s, ": ", Fehler_mittelwert(x))


# statische Methode
# Daten aus statischen generieren
T1s, T2s, T3s, T4s, T5s, T6s, T7s, T8s, ts = np.genfromtxt(
    "./content/v204_statisch.txt", unpack=True
)


# Plot mit T1 und T4
fig, ax1 = plt.subplots(label="statische Methode, T1;T4, Messing")
ax1.plot(ts, T1s, ".", label=r"$T_1$, Messing (breit)")
ax1.plot(ts, T4s, ".", label=r"$T_4$, Messing (schmal)")
ax1.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax1.legend()
fig.savefig("./build/statisch_T1_T4.pdf")

# Plot mit T5 und T8
fig, ax2 = plt.subplots(label="statische Methode, T5;T8")
ax2.plot(ts, T5s, ".", label=r"$T_5$, Aluminium")
ax2.plot(ts, T8s, ".", label=r"$T_8$, Edelstahl")
ax2.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax2.legend()
fig.savefig("./build/statisch_T5_T8.pdf")


# Plot Temperaturdifferenzen T2-T1, T7-T8
fig, (ax1, ax2) = plt.subplots(
    2, 1, label="Temperaturdifferenz T2-T1, T7-T8", layout="constrained"
)
ax1.plot(ts, T2s - T1s, ".", label=r"$\increment T_{2,1} = T_2-T_1$, Messing (breit)")
ax2.plot(ts, T7s - T8s, ".", label=r"$\increment T_{7,8} = T_7-T_8$, Edelstahl")
ax1.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$\increment T / \mathrm{°C}$",
)
ax1.legend()

######### nach Testat eingefügt #################################################################################
fig_post, (ax_post, ax_postdiff) = plt.subplots(2,1, label="statische Methode, T1;T2, Messing", layout = "constrained")
ax_post.plot(ts, T1s, ".", label=r"$T_1$, Messing (breit)")
ax_post.plot(ts, T2s, ".", label=r"$T_2$, Messing (breit)")
ax_post.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax_post.legend()

ax_postdiff.plot(ts, T2s - T1s, ".", label=r"$\increment T_{2,1} = T_2-T_1$, Messing (breit)")
ax_postdiff.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$\increment T / \mathrm{°C}$",
)
ax_postdiff.legend()

fig_post.savefig("./build/statisch_T1_T2.png")
#################################################################################################################

ax2.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$\increment T / \mathrm{°C}$",
)
ax2.legend()
fig.savefig("./build/Temperaturdifferenz.pdf")


# Auswertung
# Wärmestrom Messing (breit) berechnen
print("Messing breit t = 100s: ", Waermestrom(120, (0.012 * 0.004), 5.13))
print("Messing breit t = 200s: ", Waermestrom(120, (0.012 * 0.004), 3.71))
print("Messing breit t = 400s: ", Waermestrom(120, (0.012 * 0.004), 2.67))
print("Messing breit t = 600s: ", Waermestrom(120, (0.012 * 0.004), 2.43))
print("Messing breit t = 750s: ", Waermestrom(120, (0.012 * 0.004), 2.39))

# Wärmestrom Messing (schmal)
print("Messing schmal t = 100s: ", Waermestrom(120, (0.07 * 0.004), 5.6))
print("Messing schmal t = 200s: ", Waermestrom(120, (0.07 * 0.004), 4.0))
print("Messing schmal t = 400s: ", Waermestrom(120, (0.07 * 0.004), 3.1))
print("Messing schmal t = 600s: ", Waermestrom(120, (0.07 * 0.004), 2.94))
print("Messing schmal t = 750s: ", Waermestrom(120, (0.07 * 0.004), 2.9))

# Wärmestrom Aluminium
print("Aluminium t = 100s: ", Waermestrom(237, (0.012 * 0.004), 3.09))
print("Aluminium t = 200s: ", Waermestrom(237, (0.012 * 0.004), 1.85))
print("Aluminium t = 400s: ", Waermestrom(237, (0.012 * 0.004), 1.33))
print("Aluminium t = 600s: ", Waermestrom(237, (0.012 * 0.004), 1.26))
print("Aluminium t = 750s: ", Waermestrom(237, (0.012 * 0.004), 1.23))

# Wärmestrom Edelstahl
print("Edelstahl t = 100s: ", Waermestrom(15, (0.012 * 0.004), 9.37))
print("Edelstahl t = 200s: ", Waermestrom(15, (0.012 * 0.004), 10.0))
print("Edelstahl t = 400s: ", Waermestrom(15, (0.012 * 0.004), 8.86))
print("Edelstahl t = 600s: ", Waermestrom(15, (0.012 * 0.004), 8.26))
print("Edelstahl t = 750s: ", Waermestrom(15, (0.012 * 0.004), 8.06))


#
# dynamisch 80s
# Daten aus dynamischen Versuch, 80s generieren
T1d, T2d, _, _, T5d, T6d, _, _, td = np.genfromtxt(
    "./content/dynamisch80s.txt", unpack=True
)


# dynamisch 80s plotten
# Messing 80s plotten
fig, ax = plt.subplots(label="dynamische Methode 80s")
ax.plot(td, T1d, ".", label="dynamisch T1, Messing")
ax.plot(td, T2d, ".", label="dynamisch T2, Messing")
ax.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax.legend()
fig.savefig("./build/dynamisch_T1_T2.pdf")

#######nachträglich eingefügt ################
fig.savefig("./build/dynamisch_T1_T2.png")
##############################################

# Aluminium 80s plotten
fig, ax = plt.subplots()
ax.plot(td, T5d, ".", label="T5 dynamisch 80s")
ax.plot(td, T6d, ".", label="T6 dynamisch 80s")
ax.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax.legend()
fig.savefig("./build/dynamisch_T5_T6.pdf")


# peaks Messing finden
max_peaks_T1d, _ = find_peaks(T1d, distance=10)
min_peaks_T1d, _ = find_peaks(-T1d, distance=10)
#
max_peaks_T2d, _ = find_peaks(T2d, distance=10)
min_peaks_T2d, _ = find_peaks(-T2d, distance=10)
# Minima und Maxima zur veranschaulichung alternativ plotten
# ax.plot(td[max_peaks_T1d], T1d[max_peaks_T1d], "kx")
# ax.plot(td[min_peaks_T1d], T1d[min_peaks_T1d], "bx")
# ax.plot(td[max_peaks_T2d], T2d[max_peaks_T2d], "rx")
# ax.plot(td[min_peaks_T2d], T2d[min_peaks_T2d], "gx")
# fig.savefig("build/dynamische_80s.pdf")

# Amplituden Messing berechnen
A1 = np.zeros(9)
for i in range(9):
    A1[i] = 0.5 * (
        T1d[max_peaks_T1d[i + 1]]
        - (T1d[min_peaks_T1d[i + 1]] - T1d[min_peaks_T1d[i]]) * 0.5
        - T1d[min_peaks_T1d[i]]
    )
#
A2 = np.zeros(9)
for i in range(9):
    A2[i] = 0.5 * (
        T2d[max_peaks_T2d[i + 1]]
        - (T2d[min_peaks_T2d[i + 1]] - T2d[min_peaks_T2d[i]]) * 0.5
        - T2d[min_peaks_T2d[i]]
    )


# Phasendifferenz Messing berechnen
dt12 = np.zeros(21)
for i in range(10):
    dt12[i] = td[min_peaks_T1d[i]] - td[min_peaks_T2d[i]]
for i in range(10, 21):
    dt12[i] = td[max_peaks_T1d[i - 10]] - td[max_peaks_T2d[i - 10]]


# Peaks Aluminium finden
max_peaks_T5d, _ = find_peaks(T5d, distance=10)
min_peaks_T5d, _ = find_peaks(-T5d, distance=10)
#
max_peaks_T6d, _ = find_peaks(T6d, distance=10)
min_peaks_T6d, _ = find_peaks(-T6d, distance=10)

# Amplituden Aluminium berechnen
A5 = np.zeros(9)
for i in range(9):
    A5[i] = 0.5 * (
        T5d[max_peaks_T5d[i + 1]]
        - (T5d[min_peaks_T5d[i + 1]] - T5d[min_peaks_T5d[i]]) * 0.5
        - T5d[min_peaks_T5d[i]]
    )
#
A6 = np.zeros(9)
for i in range(9):
    A6[i] = 0.5 * (
        T6d[max_peaks_T6d[i + 1]]
        - (T6d[min_peaks_T6d[i + 1]] - T6d[min_peaks_T6d[i]]) * 0.5
        - T6d[min_peaks_T6d[i]]
    )

# Phasendifferenz Aluminium berechnen
dt56 = np.zeros(21)
for i in range(10):
    dt56[i] = td[min_peaks_T5d[i]] - td[min_peaks_T6d[i]]
for i in range(10, 21):
    dt56[i] = td[max_peaks_T5d[i - 10]] - td[max_peaks_T6d[i - 10]]


# Auswertung
# Ausgabe der Ergebnisse

# Messing
print("Kappa für Messing: ", kappa(8520, 385, A2, A1, dt12))
Wert("Amplitude Messing fern", A1)
Wert("Amplitude Messing nah", A2)
Wert("Phasendifferenz Messing (80s)", dt12)

# Aluminium
print("Kappa für Aluminium: ", kappa(2800, 830, A6, A5, dt56))
Wert("Amplitude Aluminium fern", A5)
Wert("Amplitude Aluminium nah", A6)
Wert("Phasendifferenz Aluminium (80s)", dt56)


# dynamisch 200s
# Daten dynamisch 200s
_, _, _, _, _, _, T7d, T8d, td = np.genfromtxt(
    "./content/dynamisch200s.txt", unpack=True
)

# dynamisch 200s plotten
fig, ax = plt.subplots(label="dynamische Methode 200s")
ax.plot(td, T7d, ".", label=r"$T_7$, Edelstahl")
ax.plot(td, T8d, ".", label=r"$T_8$, Edelstahl")
ax.set(
    xlabel=r"$t / \mathrm{s}$",
    ylabel=r"$T / \mathrm{°C}$",
)
ax.legend()
fig.savefig("./build/dynamisch_T7_T8.pdf")

# Edelstahl
# Peaks Edelstahl finden
max_peaks_T8d, _ = find_peaks(T8d, distance=10)
min_peaks_T8d, _ = find_peaks(-T8d, distance=10)
#
max_peaks_T7d, _ = find_peaks(T7d, distance=10)
min_peaks_T7d, _ = find_peaks(-T7d, distance=10)

# Amplituden Edelstahl berechnen
A8 = np.zeros(4)
for i in range(4):
    A8[i] = 0.5 * (
        T8d[max_peaks_T8d[i + 1]]
        - (T8d[min_peaks_T8d[i + 1]] - T8d[min_peaks_T8d[i]]) * 0.5
        - T8d[min_peaks_T8d[i]]
    )
#
A7 = np.zeros(4)
for i in range(4):
    A7[i] = 0.5 * (
        T7d[max_peaks_T7d[i + 1]]
        - (T7d[min_peaks_T7d[i + 1]] - T7d[min_peaks_T7d[i]]) * 0.5
        - T7d[min_peaks_T7d[i]]
    )

# Phasendifferenz Edelstahl berechnen
dt78 = np.zeros(9)
for i in range(4):
    dt78[i] = td[min_peaks_T7d[i]] - td[min_peaks_T8d[i]]
for i in range(4, 9):
    dt78[i] = td[max_peaks_T8d[i - 4]] - td[max_peaks_T7d[i - 4]]

# Edelstahl
print("Kappa für Edelstahl: ", kappa(8000, 400, A7, A8, dt78))
Wert("Amplitude Edelstahl fern", A8)
Wert("Amplitude Edelstahl nah", A7)
Wert("Phasendifferenz Edeltstahl (200s)", dt78)
