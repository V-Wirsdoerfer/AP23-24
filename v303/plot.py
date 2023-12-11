from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from scipy.stats import sem
from scipy.optimize import curve_fit

# Daten generieren
phi_d_n, A_n = np.genfromtxt("./content/no_noise.txt", unpack=True)
phi_d_m, A_m = np.genfromtxt("./content/mit_noise.txt", unpack=True)
x_cm, U = np.genfromtxt("./content/photodetektor.txt", unpack=True)

phi_n = phi_d_n * np.pi / 180
phi_m = phi_d_m * np.pi / 180
x = x_cm / 100


# Funktionen für curve_fit
def cosinus(x, a, b, c):
    return a * np.cos(x - b) + c

def photo(x, a, b):
    return (a / x**2) + b


# Rohdaten plotten (no noise)
fig1, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(phi_n, A_n, "x", label="Amplituden ohne Rauschen")

# curve_fit (no noise)
params_n, cov_n = curve_fit(
    cosinus,
    phi_n,
    A_n,
)
x = np.linspace(0, 6.6)
ax1.plot(
    x,
    cosinus(
        x,
        *params_n,
    ),
    "r",
    label="curve fit",
)
ax1.set(xlabel=r"$\phi / rad$", ylabel=r"$A / V$", yscale="linear")
ax1.legend()


# Rohdaten plotten (mit noise)
fig2, ax2 = plt.subplots(1, 1, layout="constrained")
ax2.plot(phi_m, A_m, "x", label="Amplituden mit Rauschen")

# curve_fit (mit noise)
params_m, cov_m = curve_fit(cosinus, phi_m, A_m)
ax2.plot(x, cosinus(x, *params_m), "r", label="curve fit")
ax2.set(
    xlabel=r"$\phi / rad$",
    ylabel=r"$A / V$"
    # yscale = "linear"
)
ax2.legend()


# Rohdaten plotten (Photodetektor)
fig3, ax3 = plt.subplots(1, 1, layout="constrained")
ax3.plot(x_cm, U, "x", label="Photodetektor")

# curve_fit (Photodetektor)
params_p, cov_p = curve_fit(photo, x_cm, U)
x = np.linspace(0.1, 90)
ax3.plot(x, photo(x, *params_p), "r", label="curve fit")
ax3.set(
    xlabel=r"$x / m$",
    ylabel=r"$U / V$",
    # xlim = (0, 90),
    ylim=(-0.5, 8.5),
)
ax3.legend()


# save
fig1.savefig("build/no_noise.pdf")
fig2.savefig("build/mit_noise.pdf")
fig3.savefig("build/photodetektor.pdf")


#Ausgabe
def error(cov):
    return np.sqrt(np.diag(cov))

er_n= error(cov_n)[0]
er_m = error(cov_m)[0]
U_n = ufloat(params_n[0], er_n)*np.pi/2 
U_m= ufloat(params_m[0], er_m)*np.pi/2
print("parameter vom Cosinus fit: \nohne Noise:\t2/pi U_0: ", U_n,  
      "\nmit Noise:\t2/pi U_0: ", U_m)


# plt.show()
