# from cProfile import label
# from distutils.command import sdist
# from turtle import distance
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


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
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax1.legend()
fig.savefig("./build/statisch_T1_T4.pdf")


# Plot mit T5 und T8
fig, ax2 = plt.subplots(label="statische Methode, T5;T8")
ax2.plot(ts, T5s, ".", label=r"$T_5$, Aluminium")
ax2.plot(ts, T8s, ".", label=r"$T_8$, Edelstahl")
ax2.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax2.legend()
fig.savefig("./build/statisch_T5_T8.pdf")


# Plot Temperaturdifferenzen T2-T1, T7-T8
fig, (ax1, ax2) = plt.subplots(
    2, 1, label="Temperaturdifferenz T2-T1, T7-T8", layout="constrained"
)
ax1.plot(ts, T2s - T1s, ".", label=r"$\Delta T_{2,1} = T_2-T_1$, Edelstahl")
ax2.plot(ts, T7s - T8s, ".", label=r"$\Delta T_{7,8} = T_7-T_8$, Edelstahl")
ax1.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax1.legend()
ax2.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax2.legend()
fig.savefig("./build/Temperaturdifferenz.pdf")


#
# dynamisch 80s
# Daten aus dynamischen Versuch, 80s generieren
T1d, T2d, _, _, T5d, T6d, _, _, td = np.genfromtxt(
    "./content/dynamisch80s.txt", unpack=True
)


# dynamisch 80s plotten
fig, ax = plt.subplots(label="dynamische Methode 80s")
ax.plot(td, T1d, ".", label="dynamisch T1, Messing")
ax.plot(td, T2d, ".", label="dynamisch T2, Messing")
ax.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax.legend()
fig.savefig("./build/dynamisch_T1_T2.pdf")

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
Wert("Amplitude Messing fern", A1)
Wert("Amplitude Messing nah", A2)
Wert("Phasendifferenz Messing (80s)", dt12)

# Aluminium
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
    xlabel=r"$t$",
    ylabel=r"$T$",
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
dt78 = np.zeros(21)
for i in range(4):
    dt78[i] = td[min_peaks_T7d[i]] - td[min_peaks_T8d[i]]
for i in range(4, 9):
    dt78[i] = td[max_peaks_T7d[i - 4]] - td[max_peaks_T8d[i - 4]]

# Edelstahl
Wert("Amplitude Edelstahl fern", A8)
Wert("Amplitude Edelstahl nah", A7)
Wert("Phasendifferenz Edeltstahl (200s)", dt78)