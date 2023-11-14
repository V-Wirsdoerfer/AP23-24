import matplotlib.pyplot as plt
import numpy as np

# Daten aus statischen generieren
T1s, T2s, T3s, T4s, T5s, T6s, T7s, T8s, ts = np.genfromtxt(
    "./content/v204_statisch.txt", unpack=True
)

# Plot mit T1 und T4
fig, ax1 = plt.subplots(label="statische Methode, T1;T4")

ax1.plot(ts, T1s, label="T1")
ax1.plot(ts, T4s, label="T4")

ax1.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax1.legend()

# Plot mit T5 und T8
fig, ax2 = plt.subplots(label="statische Methode, T5;T8")

ax2.plot(ts, T5s, label="T5")
ax2.plot(ts, T8s, label="T8")

ax2.set(
    xlabel=r"t",
    ylabel=r"T",
)
ax2.legend()


# Daten aus dynamischen Versuch, 80s generieren
T1d, T2d, _, _, _, _, _, _, td = np.genfromtxt(
    "./content/dynamisch80s.txt", unpack=True
)

# dynamisch 80s plotten
fig, ax = plt.subplots(label="dynamische Methode 80s")

ax.plot(td, T1d, label="dynamisch T1, Messing")
ax.plot(td, T2d, label="dynamisch T2, Messing")
ax.legend()


# Daten dynamisch 100s
_, _, _, _, _, _, T7d, T8d, td = np.genfromtxt("./content/dynamisch200s.txt", unpack=True)

#dynamisch 200s plotten
fig, ax = plt.subplots(label="dynamische Methode 200s")

ax.plot(td, T7d, label = "dynamisch T7, Edelstahl")
ax.plot(td, T8d, label="dynamisch T8, Edelstahl")
ax.legend()

plt.show()
# fig.savefig("build/plot.pdf")
