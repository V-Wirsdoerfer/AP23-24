from pyexpat.errors import XML_ERROR_FINISHED
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

# Daten generieren
x_Kreis_e, Dx_Kreis_e = np.genfromtxt("./content/Kreis_einseitig.txt", unpack=True)
x_Kreis_b_n, Dx_Kreis_b_n = np.genfromtxt(
    "./content/Kreis_beidseitig_nah.txt", unpack=True
)
x_Kreis_b_f, Dx_Kreis_b_f = np.genfromtxt(
    "./content/Kreis_beidseitig_fern.txt", unpack=True
)
x_Quadrat_e, Dx_Quadrat_e = np.genfromtxt(
    "./content/Quadrat_einseitig.txt", unpack=True
)
x_Quadrat_b_n, Dx_Quadrat_b_n = np.genfromtxt(
    "./content/Quadrat_beidseitig_nah.txt", unpack=True
)
x_Quadrat_b_f, Dx_Quadrat_b_f = np.genfromtxt(
    "./content/Quadrat_beidseitig_fern.txt", unpack=True
)


# Messgrößen definieren
F_einseitig = ufloat(9.81 * 0.5004, 9.81 * 0.00005)
F_beidseitig = ufloat(9.81 * 0.9993, 9.81 * 0.00005)
L_K_e = ufloat(0.535, 0.0005)
L_K_b = ufloat(0.551, 0.0005)
L_Q_e = ufloat(0.532, 0.0005)
L_Q_b = ufloat(0.555, 0.0005)
I_K = ufloat(np.pi * (0.01**4) / 64, np.pi * (0.000005**4) / 64)
I_Q = ufloat((0.01**4) / 12, (0.000005**4) / 12)


# Funktionen zum rechnen
def Dx_beidseitig_nah(x, E, L, F, I):
    return (F / 48 * E * I) * (3 * (L**2) * x - 4 * x**3)


def Dx_beidseitig_fern(x, E, L, F, I):
    return (F / 48 * E * I) * (4 * x**3 - 12 * L * x**2 + 9 * (L**2) * x - L**3)


def Dx_einseitig(x, E, L, F, I):
    return (F / 2 * E * I) * (L * x**2 - (x**3) / 3)


# linearisieren
# einseitig linearisieren
lin_x_K_e = L_K_e.nominal_value * x_Kreis_e**2 - (x_Kreis_e**3) / 3
lin_x_Q_e = L_Q_e.nominal_value * x_Quadrat_e**2 - (x_Quadrat_e**3) / 3


# beidseitig linearisieren
lin_x_K_b_n = 3 * (L_K_b.nominal_value**2) * x_Kreis_b_n - 4 * (x_Kreis_b_n**3)
lin_x_K_b_f = (
    4 * (x_Kreis_b_f**3)
    - 12 * L_K_b.nominal_value * (x_Kreis_b_f**2)
    + 9 * (L_K_b.nominal_value**2) * x_Kreis_b_f
    - L_K_b.nominal_value**3
)

lin_x_Q_b_n = 3 * (L_Q_b.nominal_value**2) * x_Quadrat_b_n - 4 * (x_Quadrat_b_n**3)
lin_x_Q_b_f = (
    4 * (x_Quadrat_b_f**3)
    - 12 * L_Q_b.nominal_value * (x_Quadrat_b_f**2)
    + 9 * (L_Q_b.nominal_value**2) * x_Quadrat_b_f
    - L_Q_b.nominal_value**3
)


# plots erstellen
# kreisförmig einseitig
fig1, ax1 = plt.subplots()
ax1.plot(lin_x_K_e, Dx_Kreis_e, "rx", label="kreisförmig, einseitig")
ax1.set(xlabel=r"$Lx^2 - \frac{x^3}{3}$")
ax1.legend()
fig1.savefig("./build/kreis_e.pdf")

# kreisförmig beidseitig nah
fig2, ax2 = plt.subplots()
ax2.plot(
    lin_x_K_b_n,
    Dx_Kreis_b_n,
    "gx",
    label=r"kreisförmig, beidseitig: $0 < x < \frac{1}{2} L$",
)
ax2.set(xlabel=r"$3L^2x - 4x^3$", ylabel=r"$D(x)$")
ax2.legend()
fig2.savefig("./build/kreis_b_n.pdf")

# kreisförmig beidseitig fern
fig3, ax3 = plt.subplots()
ax3.plot(
    lin_x_K_b_f,
    Dx_Kreis_b_f,
    "rx",
    label=r"kreisförmig, beidseitig: $\frac{1}{2} L < x < L$",
)
ax3.set(
    xlabel=r"$4 x ^3 -12 L x ^2 + 9L^2 x - L^3$",
    ylabel=r"$D(x)$",
)
ax3.legend()
fig3.savefig("./build/kreis_b_f.pdf")


# quadratisch einseitig
fig4, ax4 = plt.subplots()
ax4.plot(lin_x_Q_e, Dx_Quadrat_e, "rx", label="quadratisch, einseitig")
ax4.set(xlabel=r"$3L^2x - 4x^3$", ylabel=r"$D(x)$")
ax4.legend()

# quadratisch beidseitig nah
fig5, ax5 = plt.subplots()
ax5.plot(
    lin_x_Q_b_n,
    Dx_Quadrat_b_n,
    "bx",
    label=r"quadratisch, beidseitig: $0 < x < \frac{1}{2} L$",
)
ax5.plot(
    lin_x_Q_b_f,
    Dx_Quadrat_b_f,
    "rx",
    label=r"quadratisch, beidseitig $\frac{1}{2} L < x < L$",
)
ax5.legend()


plt.show()
