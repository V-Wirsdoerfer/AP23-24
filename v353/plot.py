from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

# import uncertainties.unumpy as unp
from scipy.optimize import curve_fit


# Daten generieren
t, U = np.genfromtxt("./content/RC-daten.txt", unpack=True)
f, A, a = np.genfromtxt("content/f-A-a.txt", unpack=True)
b = 1 / f
phi = a / b * 2 * np.pi

# Rohdaten plotten
fig1, ax1 = plt.subplots()
ax1.plot(t, U, "x", label="U(t)/t, rohdaten")
ax1.set(
    xlabel=r"t/ms",
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
uncertainties = np.sqrt(np.diag(covariance_matrix))
print("Steigung RC =", params[0], "\u00b1", uncertainties[0])

ax1.legend()


# b)
# halblogarithmischer Plot
omega = f * 2 * np.pi


def A_w(omega, a, b):
    return b / np.sqrt(1 + omega**2 * a**2)


fig2, ax2 = plt.subplots()
ax2.plot(omega, A, "x", label=r"$A(v_i),(v_i)$")
ax2.set(
    xlabel=r"$v_i/s^{-1}$",
    ylabel=r"$A(v_i)/mV$",
    yscale="log",
    ylim=(0, 1),
)

a = -15
b = 5
b = 40000000
x = np.linspace(10, 650000, 1000)
params, covariance_matrix = curve_fit(A_w, omega, A, p0=(a, b))
ax2.plot(x, A_w(x, *params), label="curve fit")
ax2.plot(x, (b / np.sqrt(1 + x**2 * a**2)), label="per Hand")

ax2.legend()
print(params)

fig3, ax3 = plt.subplots()
ax3.plot(omega, phi, "x", label="Phasenverschiebung")
ax3.legend()
ax3.set(
    xlabel=r"$\nu_i$",
    ylabel=r"$\phi(\nu_i)$",
    yscale="log",
)

x = 2 * np.pi * np.linspace(30, 100000,)
a = -0.0000000052
# b=5e-0
# c=-1e4

a, b, c = 1, 1, 1


def phi_w(omega, a, b, c):
    return np.arctan(-omega * a) * b + c


params, covariance_matrix = curve_fit(phi_w, omega, a, p0=(a, b, c))
ax3.plot(
    x,
    np.cos(-x * params[0]) / np.sin(-x * params[0]) * params[1] + params[2],
    label="per hand",
)
#ax3.plot(x, np.cos(-x * a) / np.sin(-x * a) * b + c, label="per hand")
print("params ax3: ", *params)
ax3.legend()
plt.show()
