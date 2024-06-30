from turtle import color
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from uncertainties import unumpy as unp


### Naturkonstanten
h = 4.1336e-15  # in eVs
c = 3e8         # in m/s
alpha = 7.297e-3    # Feinstrukturkonstanten, einheitenlos
d = 201.4e-12   #in m

### Daten generieren

theta_Bragg, Rate_Bragg = np.genfromtxt("BraggBed.txt", unpack=True)
theta_Cu, Rate_Cu = np.genfromtxt("EmissionCu.txt", unpack=True)
theta_Zn, Rate_Zn = np.genfromtxt("EmissionZn.txt", unpack=True)
theta_Ga, Rate_Ga = np.genfromtxt("EmissionGa.txt", unpack=True)
theta_Br, Rate_Br = np.genfromtxt("EmissionBr.txt", unpack=True)
theta_Sr, Rate_Sr = np.genfromtxt("EmissionSr.txt", unpack=True)
theta_Zr, Rate_Zr = np.genfromtxt("EmissionZr.txt", unpack=True)


# 2 theta in theta umrechnen
theta_Cu /= 2
theta_Zn /= 2
theta_Ga /= 2
theta_Br /= 2
theta_Sr /= 2
theta_Zr /= 2


# print(theta_Bragg, Rate_Bragg)


### Funktionen definieren
def Wellenlänge(theta, n=1):
    """Wird mithilfe der Braggschen Bedingung berechnet."""
    theta = theta / 360 * (2 * np.pi)
    return 2 / n * d * np.sin(theta)


def Steigung(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[0], err[0])


def Achsenabschnitt(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[1], err[1])


def half_height(ax, Steigbereich, Data, name):
    """
    Bestimmt die Stelle der halben Höhe und zeichnet diese ein
    returns Stelle der Halben Höhe
    """
    params, cov = np.polyfit(Steigbereich, Data, deg=1, cov=True)
    x = np.linspace(Steigbereich[0], Steigbereich[-1])
    ax.plot(
        x, Steigung(params, cov).n * x + Achsenabschnitt(params, cov).n, label=f"{name}"
    )
    Mitte = (min(x) + max(x)) / 2
    Half_height = Steigung(params, cov) * Mitte + Achsenabschnitt(params, cov)
    ax.errorbar(
        Mitte,
        Steigung(params, cov).n * Mitte + Achsenabschnitt(params, cov).n,
        #    yerr=unp.std_devs(Half_height),
        fmt="x",
        label="halbe Höhe",
    )

    print(f"Die Halbe Höhe von {name} befindet sich bei ({Mitte}, {Half_height})")
    return Mitte


def half_height_Cu(No1, No2, No3, No4, max):
    """
    Nox beschreibt den Index im Array
    max ist der Maximalwert der Funktion
    Returns the x-Value of the intsect
    """
    m_l = (Rate_Cu[No2] - Rate_Cu[No1]) / (theta_Cu[No2] - theta_Cu[No1])
    b_l = Rate_Cu[No1] - m_l * theta_Cu[No1]
    m_r = (Rate_Cu[No4] - Rate_Cu[No3]) / (theta_Cu[No4] - theta_Cu[No3])
    b_r = Rate_Cu[No3] - m_r * theta_Cu[No3]

    intsect_l = (max / 2 - b_l) / m_l
    intsect_r = (max / 2 - b_r) / m_r
    return [intsect_l, intsect_r]


def E_K(theta, n=1):
    """
    Berechnet die Energie zu einem festen Winkel
    alternative Ordnung n wählbar
    """
    return h*c/Wellenlänge(theta, n)


def sigma_K(Z, E_K):
    """
    Z ist Ordnungszahl und E_K ist die Energie der K-Kante"""
    return Z - unp.sqrt(E_K / Ryd - alpha**2 * Z**4 / 4)


### Daten plotten

# Bragbedingung
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(theta_Bragg, Rate_Bragg, ".", label="Braggbedingung")
print(
    "Maximaler Bragg Winkel:",
    theta_Bragg[Rate_Bragg == max(Rate_Bragg)],
    "mit Zählrate von: ",
    max(Rate_Bragg),
)
ax1.plot(27.7, max(Rate_Bragg), "o", label="Braggwinkel")
ax1.vlines(x=27.7, ymin=0, ymax=max(Rate_Bragg), color="orange", linestyle="dotted")

ax1.set(
    ylim=[0, 89],
    xlabel=r"$2 \cdot \theta/°$",
    ylabel="Zählrate",
)
ax1.legend()
fig1.savefig("build/BraggBed.pdf")

# Emissionsspektrum Cu-Röntgenröhre
fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(theta_Cu, Rate_Cu, ".", label="Emissionsspektrum Cu")
# K-Linien
ax2.vlines(
    theta_Cu[81], 0, Rate_Cu[81], linestyle="dotted", color="orange", label=r"$K_\beta$"
)
ax2.vlines(
    theta_Cu[92], 0, Rate_Cu[92], linestyle="dotted", color="green", label=r"$K_\alpha$"
)

print(
    f"Winkel der K beta Linie: {theta_Cu[81]}°",
    f"\nWinkel der K alpha Linie: {theta_Cu[92]}°",
)
E_K_beta = E_K(theta_Cu[81])
print(f"Energie K beta = {E_K_beta}eV")
E_K_alpha = E_K(theta_Cu[92])
print(f"Energie K alpha = {E_K_alpha}eV")

# Bremsspektrum
ax2.fill_between(
    theta_Cu[0:79],
    Rate_Cu[0:79] - 60,
    Rate_Cu[0:79] + 60,
    alpha=0.2,
    color="b",
    label="Bremsberg",
)
ax2.fill_between(
    theta_Cu[83:90], Rate_Cu[83:90] - 60, Rate_Cu[83:90] + 60, alpha=0.2, color="b"
)
ax2.fill_between(
    theta_Cu[95:], Rate_Cu[95:] - 60, Rate_Cu[95:] + 60, alpha=0.2, color="b"
)

ax2.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax2.legend()
fig2.savefig("build/EmissionCu.pdf")

# zoom
ax2.set(xlim=[18.5, 24])
fig2.savefig("build/EmissionCU_zoom.pdf")

# Halfheigt
ax2.plot(theta_Cu[78:84], Rate_Cu[78:84], "b")
ax2.plot(theta_Cu[89:96], Rate_Cu[89:96], "b")
intsect_beta = half_height_Cu(79, 80, 81, 82, Rate_Cu[81])
intsect_alpha = half_height_Cu(90, 91, 93, 94, Rate_Cu[92])

ax2.hlines(
    y=Rate_Cu[81] / 2,
    xmin=intsect_beta[0],
    xmax=intsect_beta[1],
    color="r",
    linestyle="dotted",
    label=r"Halbe Höhe von $K_\beta$",
)
ax2.hlines(
    y=Rate_Cu[92] / 2,
    xmin=intsect_alpha[0],
    xmax=intsect_alpha[1],
    color="k",
    linestyle="dotted",
    label=r"Halbe Höhe von $K_beta$",
)
ax2.legend()
fig2.savefig("build/Halfheight.pdf")
print(f"Breite bei halber Höhe beta ={intsect_beta[1]-intsect_beta[0]}°")
print(f"Breite bei halber Höhe alpha={intsect_alpha[1]-intsect_alpha[0]}°")

# Grenzwinkel/ max E
ax2.set(xlim=[3.5, 11], ylim=[-50, 200])
fig2.savefig("build/Grenzwinkel.pdf")

print(f"minimale Wellenlänge = {Wellenlänge(5, 1)}m")
print(f"maximale Energie {E_K(5)} eV")

# Absorptionsspektrum Zn Probe
fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(theta_Zn, Rate_Zn, "x", label="Absorptionsspektrum Zn")
theta_K_Zn = half_height(ax3, theta_Zn[2:6], Rate_Zn[2:6], "Zink")

ax3.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax3.legend()
fig3.savefig("build/AbsorptionZn.pdf")

# Absorptionsspektrum Ga Probe
fig4, ax4 = plt.subplots(layout="constrained")
ax4.plot(theta_Ga, Rate_Ga, "x", label="Absorptionsspektrum Ga")
theta_K_Ga = half_height(ax4, theta_Ga[2:7], Rate_Ga[2:7], "Gallium")

ax4.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax4.legend()
fig4.savefig("build/AbsorptionGa.pdf")

# Absorptionsspektrum Br Probe
fig5, ax5 = plt.subplots(layout="constrained")
ax5.plot(theta_Br, Rate_Br, "x", label="Absorptionsspektrum Br")
theta_K_Br = half_height(ax5, theta_Br[2:7], Rate_Br[2:7], "Brom")

ax5.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax5.legend()
fig5.savefig("build/AbsorptionBr.pdf")

# Absorptionsspektrum Sr Probe
fig6, ax6 = plt.subplots(layout="constrained")
ax6.plot(theta_Sr, Rate_Sr, "x", label="Absorptionsspektrum Sr")
theta_K_Sr = half_height(ax6, theta_Sr[2:7], Rate_Sr[2:7], "Gallium")

ax6.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax6.legend()
fig6.savefig("build/AbsorptionSr.pdf")

# Absorptionsspektrum Zr Probe
fig7, ax7 = plt.subplots(layout="constrained")
ax7.plot(theta_Zr, Rate_Zr, "x", label="Absorptionsspektrum Zr")
theta_K_Zr = half_height(ax7, theta_Zr[3:7], Rate_Zr[3:7], "Zirkonium")

ax7.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax7.legend()
fig7.savefig("build/AbsorptionZr.pdf")


### Abschirmkonstanten Kupfer
E_K_abs = 8.997e3  # in eV
Z_Kupfer = 29  # Ordnungszahl
Ryd = 13.6  # Rydbergenergie
m = 2
l = 3

sigma1 = Z_Kupfer - np.sqrt(E_K_abs / Ryd)
print(f"sigma1 = {sigma1}")

sigma2 = -m * np.sqrt((-E_K_alpha + E_K_abs) / Ryd) + Z_Kupfer
print(f"sigma2 = {sigma2}")

sigma3 = -l * np.sqrt((-E_K_beta + E_K_abs) / Ryd) + Z_Kupfer
print(f"sigma3 = {sigma3}")


### Energien der K-Kanten für die fünf Metalle
theta_K = np.asarray([theta_K_Br, theta_K_Ga, theta_K_Sr, theta_K_Zn, theta_K_Zr])
E_K_Metalle = E_K(theta_K)
name = ["Brom", "Gallium", "Strontium", "Zink", "Zirkonium"]
Z_Metalle = np.asarray([35, 31, 38, 30, 40])
sigma_K_Metalle = sigma_K(Z_Metalle, E_K_Metalle)
for i in range(5):
    print(f"Theta der halben Höhe von {name[i]} ist {theta_K[i]}°")
    print(f"Absorptionsenergie der K-Kante von {name[i]} beträgt {E_K_Metalle[i]} eV.")
    print(f"Die Abschirmkonstante von {name[i]} beträgt {sigma_K_Metalle[i]}")


### Moseleysches-Gesetz
fig8, ax8 = plt.subplots(layout="constrained")
ax8.plot(Z_Metalle, np.sqrt(E_K_Metalle), ".", label="Absorptionsenergien der K-Kante")
ax8.set(
    xlabel="Z",
    ylabel = r"$\sqrt{E_K}/\sqrt{eV}$",
)
params8, cov8 = np.polyfit(Z_Metalle, np.sqrt(E_K_Metalle), deg=1, cov=True)
x8=np.linspace(min(Z_Metalle), max(Z_Metalle))
m8 = Steigung(params8, cov8)
b8 = Achsenabschnitt(params8, cov8)
ax8.plot(x8, m8.n*x8+b8.n, label="Ausgleichsgerade" )
ax8.legend()
fig8.savefig("build/Moseley.pdf")

print(f"Mooseley: Steigung = {m8}, Achsenabschnitt = {b8} \nRydberg-Energie: {m8**2}")



### Abweichung von Theoriewerten
def theta_aus_E(E_K, n=1):
    return np.arcsin(n * c * h/(2*d * E_K))







