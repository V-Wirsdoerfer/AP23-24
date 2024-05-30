import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat


#Daten generieren
Winkel, Intensität_senkrecht, drittel_senkrecht, Intensität_parallel, drittel_parallel, Fehler_senkrecht, Fehler_parallel = np.genfromtxt("Messwerte.txt", unpack=True)  #in °, A, A, A, A

I0 = ufloat(0.48e-3, 0.01e-3)   #in A


#Manche Werte durch 3 teilen
drittel_senkrecht = np.asarray(drittel_senkrecht, dtype=bool)
drittel_parallel  = np.asarray(drittel_parallel,  dtype=bool)

Intensität_senkrecht[drittel_senkrecht] /= 3
Intensität_parallel[drittel_parallel]   /= 3


print("Intensität parallel:\n", np.round(Intensität_senkrecht, 7))

fig1, ax1 = plt.subplots(layout="constrained")

ax1.errorbar(Winkel,  Intensität_senkrecht, fmt="x", label="senkrecht")
ax1.errorbar(Winkel,  Intensität_parallel, fmt="x", label="parallel")
ax1.set(
    xlabel="Winkel / °",
    ylabel="Photostrom / A",
)
ax1.legend()

fig1.savefig("build/Messdaten.pdf")


#wurzel I/I0
fig2, ax2 = plt.subplots(layout="constrained")

ax2.errorbar(Winkel,  np.sqrt(Intensität_senkrecht/unp.nominal_values(I0)), fmt="x", label="senkrecht")
ax2.errorbar(Winkel,  np.sqrt(Intensität_parallel/unp.nominal_values(I0)), fmt="x", label="parallel")
ax2.set(
    xlabel="Winkel / °",
    ylabel=r"$\sqrt{\frac{I}{I_0}}$",
)
ax2.legend()

fig2.savefig("build/I0.pdf")

