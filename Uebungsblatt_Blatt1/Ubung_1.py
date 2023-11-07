import numpy as np
import matplotlib.pyplot as plt

# Daten Generieren
x, N = np.genfromtxt("Daten.txt", unpack=True)

# Daten Plotten
fig, ax1 = plt.subplots()

ax1.plot(x, N, "x", label="linear")
ax1.set(
    title="Logarithmische Auftragung",
    xlabel="d/[cm]",
    ylabel="N/Anzahl der Gamma-Quanten",
    yscale="log",
)

plt.savefig("Auswertung_log.pdf")

fig, ax2 = plt.subplots()

ax2.plot(x, N, "x")
ax2.set(
    title="Lineare Auftragung",
    xlabel=r"$d$/[cm]",
    ylabel="N/Anzahl der Gammaquanten",
    yscale="linear",
)

plt.savefig("Auswertung_linear.pdf")