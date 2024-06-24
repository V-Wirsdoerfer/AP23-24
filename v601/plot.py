import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

T = np.array([25, 141.4, 165, 165, 174, 186])

# for i in range(len(T)):
T += 273.15


def p_sät(T):
    return 5.5e7 * np.exp(-6876 / T)


def w(T):
    return 0.0029 / p_sät(T)


print("Temperaturen:\n", T, "K")
print("mittlere Weglänge w:\n", w(T), "cm")
print("Dampfdruck:\n", p_sät(T), "mbar")
print("Verhältnis a/w:", 1 / w(T))


# Diagramme sind zu zeichnen blabla sinnlos

U = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
I_blau = [
    0.00011,
    -1e-05,
    0.00045,
    0.00022,
    0.00055,
    0.00078,
    0.0019,
    0.0097,
    0.0063,
    0,
]
I_rot = [
    0.30,
    0.26,
    0.27,
    0.04,
    0.03,
    0.03,
    0.04,
    0.03,
    0.00,
    0.00,
]

fig1, (ax1, ax2) = plt.subplots(1, 2, layout="constrained")
# ax = np.ravel(ax1)
ax1.plot(U, I_blau, "bx", label="Steigung")
ax2.plot(U, I_rot, "rx", label="Steigung")
ax1.legend()
ax1.legend()
ax1.set(
    xlabel="U/V",
    ylabel=r"$\Delta I_A / mA$",
)
ax1.set(
    xlabel="U/V",
    ylabel=r"$\Delta I_A / nA$",
)
fig1.savefig("build/Steigung.pdf")


### Mittelwert Energieniveaus

delta_U = [-4.91, -4.86, -4.82, -5, -5.18]

mean_U = ufloat(np.mean(delta_U), np.std(delta_U))
print("mittlerwert Delta U = ", mean_U, "V")



### Wellenlänge berechnen
h = 4.135e-15       #in eV*s
c=3e8               #in m/s
e=1                 #in e
lamb = h * c /(e * mean_U)

print("Die Wellenlänge beträgt: ", lamb, "m" )


### Abweichung ausrechnen

delta_U = (mean_U + 4.89)/4.89
print("Die Abweichung vom Literaturwert für die Energie beträgt: ", delta_U*100)

lamb_theo = h * c /(e * 4.89)
delta_lamb = (lamb + lamb_theo)/lamb_theo
print("Die Abweichung vom Literaturwert für die Energie beträgt: ", delta_lamb*100)
print("Theoriewert lambda: ", lamb_theo, "m")