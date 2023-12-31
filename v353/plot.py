from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from scipy.stats import sem
from scipy.optimize import curve_fit

# import uncertainties.unumpy as unp

# Daten generieren
t, U_v = np.genfromtxt("./content/t-U.txt", unpack=True)
f, A, a = np.genfromtxt("content/f-A-a.txt", unpack=True)
b = 1 / f
phi = a / b * 2 * np.pi
# print("phi:", phi, "\na: ", a, "\nb: ", b)

U = U_v / 5


def uncertainties(cov):
    return np.sqrt(np.diag(covariance_matrix))


# Rohdaten plotten
fig1, ax1 = plt.subplots()
ax1.plot(t, U, "x", label=r"$U(t)$ Rohdaten")
ax1.set(
    xlabel=r"$t/s$",
    ylabel=r"$U(t)/V$",
    yscale="linear",
)

# Daten logarithmieren
U_ln = np.log(U[0:-1])

# logarithmierte Daten plotten
ax1.plot(t[0:-1], U_ln, "x", label=r"$U(t)$ logarithmiert")

# lineare regression mit logarthmierten Daten erstellen und plotten
params, covariance_matrix = np.polyfit(t[0:-1], U_ln, deg=1, cov=True)
ax1.plot(t[0:-1], t[0:-1] * params[0] + params[1], label="lineare Regresion")
m = ufloat(params[0], uncertainties(covariance_matrix)[0])
RC_a = -1 / m
print("Steigung RC =", RC_a)

ax1.legend()
fig1.savefig("build/Entladungskurve.pdf")

# b)
# halblogarithmischer Plot
# Konstanten
omega = f * 2 * np.pi

U_0 = 5


# zu fittende Funktion
def A_w(omega, a):
    return U_0 / np.sqrt(1 + omega**2 * a**2)


fig2, ax2 = plt.subplots()
ax2.plot(f, A, "x", label=r"$A(f_i),(f_i)$")

a = 0.00152
x = np.linspace(10, 140000, 1000)
params, covariance_matrix = curve_fit(A_w, f * 2 * np.pi, A, p0=(a))
ax2.plot(x, A_w(x * 2 * np.pi, *params), label="curve fit plot")
ax2.plot(
    x, (U_0 / np.sqrt(1 + (x * 2 * np.pi) ** 2 * a**2)), label="mit Werten aus a)"
)
ax2.legend()
RC_b = ufloat(params[0], uncertainties(covariance_matrix)[0])
print("RC mit Amplituden: ", RC_b)

ax2.set(
    xlabel=r"$f_i/\unit{\hertz}$",
    ylabel=r"$A(f_i)/V$",
    yscale="log",
    ylim=(0, 1),
)
fig2.savefig("build/Amplitudenspannung.pdf")


# c)

#print("c:\n\n", "phi: ", phi, "\nomega: ", omega, "\nRC: ", np.tan(phi) / -omega, "\n")

# Mittelwertfehler
# B = [2.5297e-3, 9.0056e-3, 3.373e-4, 2.53e-4, 3.376e-4]
# B_M = ufloat(np.mean(B), sem(B))
# print("Mittelwert gute Daten", B_M)
# print("Mittelwertfehler aus Guten Daten", sem(B))
# print("phi später: ", phi)


# zu fittende Funktion
def phi_w(omega, a):
    return np.arctan(-omega * a)


fig3, ax3 = plt.subplots()

params, covariance_matrix = curve_fit(phi_w, omega, a, p0=(a))

# x = np.linspace(0, 120000, 100000)
x = np.logspace(0, 5)

# ax3.plot(x, np.cos(-x * params[0]) / np.sin(-x * params[0]), label="per hand")
ax3.plot(f, phi, "x", label=r"Phasenverschiebung $\varphi$")
#ax3.plot(x, np.arctan(x), label="arctan(x)")
ax3.plot(x, -phi_w(x * 2 * np.pi, RC_a.nominal_value), label="Werte aus a)")
ax3.plot(x, -phi_w(x * 2 * np.pi, RC_b.nominal_value), label="Werte aus b)")
ax3.plot(x, -phi_w(x * 2 * np.pi, 2.5e-3), label=r"$RC=2.5e-3$")
ax3.legend()
ax3.set(
    xlabel=r"$f_i/\unit{\hertz}$",
    ylabel=r"$\varphi(f_i)$",
    xscale="log",
    # ylim=(0, 3),
    xlim=(1,1.2e+5),
    #    xlim =(-7e+5,7e+2)
)
fig3.savefig("build/phi(f).pdf")
# print("params ax3: ", *params)


# d)
# Polarplot
# phi = np.arctan(-omega*RC_a.nominal_value)
# omega = omega / (2 * np.pi) * 360
# phi = phi / (2 * np.pi) * 360


#zu fittende Funktion
def A_f():
    return -np.sin(phi) / (omega * RC_a.nominal_value)

A_f = abs(A_f())
print("phi: ", phi, "\nA_f(): ", A_f)

fig4, ax4 = plt.subplots(subplot_kw={"projection": "polar"})
ax4.plot(phi, A_f, "x")
z = np.linspace(np.pi * 0, 0.5 * np.pi, 5000)
ax4.plot(z, np.sin((z + 0.5 * np.pi) * 1))

ax4.set_rmin(0)
ax4.set_rmax(0.11)
fig4.savefig("build/polar1.pdf")

ax4.set_rmax(1)
fig4.savefig("build/polar2.pdf")


# ax4 = plt.polar(phi, A_f())
# plt.show
# def A_w(w):
#    return -np.sin(np.arctan(-w* RC_a.nominal_value))/ (w * RC_a.nominal_value)
#
# def A_curvefit(a):
#    return -np.sin(phi) / (omega * a) * U_0
#
#
#
##print("f: ", A_f(), "phi: ", phi)
#
# A_f = abs(A_f())
# w = np.linspace(0, 650000)
# ax4.plot(phi, A, "rx")
# ax4.plot(w, A_w(w), label = "theoretische A(w) Kurve")
# print("A_f(): ", A_f)
##ax4.plot(0.5 * np.pi)
## params, _ = curve_fit(A_curvefit, )
#
## ax4.plot(phi, label = "phi")
## ax4.plot(omega, U_0/np.sqrt(1+omega**2*RC_a.nominal_value **2),"r")
#
##plt.tight_layout(pad=0, h_pad=0.08, w_pad=0.08)
#
# ax4.legend()
