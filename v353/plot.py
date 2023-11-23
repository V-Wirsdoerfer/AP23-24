from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

# import uncertainties.unumpy as unp
from scipy.optimize import curve_fit


# Daten generieren
t, U_v = np.genfromtxt("./content/RC-daten.txt", unpack=True)
f, A, a = np.genfromtxt("content/f-A-a.txt", unpack=True)
b = 1 / f
phi = a / b * 2 * np.pi
# print("phi:", phi, "\na: ", a, "\nb: ", b)

U = U_v / 5


def uncertainties(cov):
    return np.sqrt(np.diag(covariance_matrix))


# Rohdaten plotten
fig1, ax1 = plt.subplots()
ax1.plot(t, U, "x", label="U(t)/t, rohdaten")
ax1.set(
    xlabel=r"t/s",
    ylabel=r"U(t)/V",
    yscale="linear",
)

# Daten logarithmieren
U_ln = np.log(U[0:-1])

# logarithmierte Daten plotten
ax1.plot(t[0:-1], U_ln, "x", label="U(t) logarithmiert")

# lineare regression mit logarthmierten Daten erstellen und plotten
params, covariance_matrix = np.polyfit(t[0:-1], U_ln, deg=1, cov=True)
ax1.plot(t[0:-1], t[0:-1] * params[0] + params[1], label="lineare Regresion")
m = ufloat(params[0], uncertainties(covariance_matrix)[0])
RC = -1 / m
print("Steigung RC =", RC)

ax1.legend()
fig1.savefig("build/Entladungskurve.pdf")

# b)
# halblogarithmischer Plot
omega = f * 2 * np.pi

U_0 = 5


def A_w(omega, a):
    return U_0 / np.sqrt(1 + omega**2 * a**2)


fig2, ax2 = plt.subplots()
ax2.plot(omega, A, "x", label=r"$A(v_i),(v_i)$")
ax2.set(
    xlabel=r"$v_i/s^{-1}$",
    ylabel=r"$A(v_i)/V$",
    yscale="log",
    ylim=(0, 1),
)

a = 0.00152
x = np.linspace(10, 650000, 1000)
params, covariance_matrix = curve_fit(A_w, omega, A, p0=(a))
ax2.plot(x, A_w(x, *params), label="curve fit plot")
ax2.plot(x, (U_0 / np.sqrt(1 + x**2 * a**2)), label="mit Werten aus a)")
ax2.legend()
RC = ufloat(params[0], uncertainties(covariance_matrix)[0])
print("RC mit Amplituden: ", RC)
fig2.savefig("build/Amplitudenspannung.pdf")


fig3, ax3 = plt.subplots()
# print("phi später: ", phi)
ax3.plot(omega, phi, "x", label="Phasenverschiebung")
ax3.legend()
ax3.set(
    xlabel=r"$\nu_i$",
    ylabel=r"$\phi(\nu_i)$",
    xscale="log",
    ylim = (0,3)
)

x = (
    2
    * np.pi
    * np.linspace(
        30,
        100000,
    )
)
a = 520
# b=5e-0
# c=-1e4

a, b, c = 1, 1, 1


def phi_w(omega, a):
    return np.arctan(-omega * a)


params, covariance_matrix = curve_fit(phi_w, omega, a, p0=(a))
ax3.plot(
    x,
    np.cos(-x * params[0]) / np.sin(-x * params[0]),
    label="per hand",
)
ax3.plot(x, np.arctan(x), label="per hand")
# print("params ax3: ", *params)
ax3.legend()
fig3.savefig("build/phi(f).pdf")
