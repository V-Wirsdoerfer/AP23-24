from cProfile import label
from turtle import distance
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


# statisch
# Daten aus statischen generieren
T1s, T2s, T3s, T4s, T5s, T6s, T7s, T8s, ts = np.genfromtxt(
    "./content/v204_statisch.txt", unpack=True
)

# Plot mit T1 und T4
fig, ax1 = plt.subplots(label="statische Methode, T1;T4, Messing")
ax1.plot(ts, T1s, label=r"T_1, Messing")
ax1.plot(ts, T4s, label=r"T_4, Messing")
ax1.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax1.legend()

# Plot mit T5 und T8
fig, ax2 = plt.subplots(label="statische Methode, T5;T8")
ax2.plot(ts, T5s, label=r"$T_5$, Aluminium")
ax2.plot(ts, T8s, label=r"$T_8$, Edelstahl")
ax2.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax2.legend()

# Temperaturdifferenzen T2-T1, T7-T8
fig, (ax1, ax2) = plt.subplots(
    2, 1, label="Temperaturdifferenz T2-T1, T7-T8", layout="constrained"
)
ax1.plot(ts, T2s - T1s, label=r"$\Delta T_{2,1} = T_2-T_1$, Edelstahl")
ax2.plot(ts, T7s - T8s, label=r"$\Delta T_{7,8} = T_7-T_8$, Edelstahl")
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


# Dynamisch 80s
# Daten aus dynamischen Versuch, 80s generieren
T1d, T2d, _, _, T5d, T6d, _, _, td = np.genfromtxt(
    "./content/dynamisch80s.txt", unpack=True
)

# dynamisch 80s plotten
fig, ax = plt.subplots(label="dynamische Methode 80s")
ax.plot(td, T1d, label="dynamisch T1, Messing")
ax.plot(td, T2d, label="dynamisch T2, Messing")


# peaks Messing finden
max_peaks_T1d, _ = find_peaks(T1d, distance=10)
ax.plot(td[max_peaks_T1d], T1d[max_peaks_T1d], "kx")
min_peaks_T1d, _ = find_peaks(-T1d, distance=10)
ax.plot(td[min_peaks_T1d], T1d[min_peaks_T1d], "bx")
#
max_peaks_T2d, _ = find_peaks(T2d, distance=10)
ax.plot(td[max_peaks_T2d], T2d[max_peaks_T2d], "rx")
min_peaks_T2d, _ = find_peaks(-T2d, distance=10)
ax.plot(td[min_peaks_T2d], T2d[min_peaks_T2d], "gx")
ax.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax.legend()
fig.savefig("build/dynamische_80s.pdf")

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
    print(T2d[max_peaks_T2d[i + 1]], T2d[min_peaks_T2d[i + 1]], T2d[min_peaks_T2d[i]])
    A2[i] = 0.5 * (
        T2d[max_peaks_T2d[i + 1]]
        - (T2d[min_peaks_T2d[i + 1]] - T2d[min_peaks_T2d[i]]) * 0.5
        - T2d[min_peaks_T2d[i]]
    )
    print(A2[i])
#print("A2: ", A2)
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
max_peaks_T6d, _ = find_peaks(T2d, distance=10)
min_peaks_T6d, _ = find_peaks(-T2d, distance=10)

# Amplituden Aluminium berechnen
A5 = np.zeros(9)
for i in range(9):
    A1[i] = 0.5 * (
        max_peaks_T5d[i + 1]
        - (min_peaks_T5d[i + 1] - min_peaks_T5d[i]) * 0.5
        - min_peaks_T5d[i]
    )
#
A6 = np.zeros(9)
for i in range(9):
    A6[i] = 0.5 * (
        max_peaks_T6d[i + 1]
        - (min_peaks_T6d[i + 1] - min_peaks_T6d[i]) * 0.5
        - min_peaks_T6d[i]
    )

# Phasendifferenz Aluminium berechnen
dt56 = np.zeros(21)
for i in range(10):
    dt56[i] = td[min_peaks_T5d[i]] - td[min_peaks_T6d[i]]
for i in range(10, 21):
    dt56[i] = td[max_peaks_T5d[i - 10]] - td[max_peaks_T6d[i - 10]]


# Fehlerrechnung
def mittelwert(x):
    return sum(x) / np.size(x)


def Fehler_mittelwert(x):
    x_M = mittelwert(x)
    sum = 0
    for i in range(np.size(x)):
        sum += (x[i] - x_M) ** 2
    x_F = np.sqrt(sum / (np.size(x) * (np.size(x) - 1)))


# Auswertung
# Mittelwert Amplitude A1
# A1 = np.array([1,2,3,4,5])
print(np.size(A1))
xA1 = mittelwert(A1)
x_FA1 = Fehler_mittelwert(A1)

print(xA1)
print(x_FA1)

#plt.show()
# a = np.arange(13)
# print(np.size(a))
# b = sum(a)
# print(a, b)

# dynamisch 200s
# Daten dynamisch 200s
_, _, _, _, _, _, T7d, T8d, td = np.genfromtxt(
    "./content/dynamisch200s.txt", unpack=True
)

# dynamisch 200s plotten
fig, ax = plt.subplots(label="dynamische Methode 200s")
ax.plot(td, T7d, label=r"$T_7$, Edelstahl")
ax.plot(td, T8d, label=r"$T_8$, Edelstahl")
ax.set(
    xlabel=r"$t$",
    ylabel=r"$T$",
)
ax.legend()

# plt.show()

# fig.savefig("build/plot.pdf")
