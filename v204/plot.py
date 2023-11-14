from cProfile import label
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
ax1.plot(ts, T1s, label="T1, Messing")
ax1.plot(ts, T4s, label="T4, Messing")
ax1.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax1.legend()

# Plot mit T5 und T8
fig, ax2 = plt.subplots(label="statische Methode, T5;T8")
ax2.plot(ts, T5s, label=r"$T5$, Aluminium")
ax2.plot(ts, T8s, label=r"$T8$, Edelstahl")
ax2.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax2.legend()

# Temperaturdifferenzen T2-T1, T7-T8
fig, (ax1, ax2) = plt.subplots(
    2, 1, label="Temperaturdifferenz T2-T1, T7-T8", layout="constrained"
)
ax1.plot(ts, T2s - T1s, label=r"$\Delta T_{st} = T2-T1$, Edelstahl")
ax2.plot(ts, T7s - T8s, label=r"$\Delta T_{st} = T7-T8$, Edelstahl")
ax1.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax1.legend()
ax2.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax2.legend()


# Dynamisch 80s
# Daten aus dynamischen Versuch, 80s generieren
T1d, T2d, _, _, _, _, _, _, td = np.genfromtxt(
    "./content/dynamisch80s.txt", unpack=True
)

# dynamisch 80s plotten
fig, ax = plt.subplots(label="dynamische Methode 80s")
ax.plot(td, T1d, label="dynamisch T1, Messing")
ax.plot(td, T2d, label="dynamisch T2, Messing")
ax.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax.legend()


# dynamisch 200s
# Daten dynamisch 200s
_, _, _, _, _, _, T7d, T8d, td = np.genfromtxt(
    "./content/dynamisch200s.txt", unpack=True
)

# dynamisch 200s plotten
fig, ax = plt.subplots(label="dynamische Methode 200s")
ax.plot(td, T7d, label=r"$T7$, Edelstahl")
ax.plot(td, T8d, label=r"$T8$, Edelstahl")
ax.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax.legend()

plt.show()

# fig.savefig("build/plot.pdf")
