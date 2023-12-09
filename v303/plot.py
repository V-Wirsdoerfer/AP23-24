from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from scipy.stats import sem
from scipy.optimize import curve_fit

#Daten generieren 
phi_d_n, A_n = np.genfromtxt("./content/no_noise.txt", unpack=True)
phi_d_m, A_m = np.genfromtxt("./content/mit_noise.txt", unpack=True)
x_cm, U = np.genfromtxt("./content/photodetektor.txt", unpack=True)

phi_n = phi_d_n * np.pi / 180
phi_m = phi_d_m * np.pi / 180
x = x_cm / 100

 
#Rohdaten plotten (no noise)
fig1, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(phi_n, A_n, "x", label="Amplituden ohne Rauschen")
ax1.set(
    xlabel = r"$\phi / rad$",
    ylabel = r"$A / V$",
    yscale = "linear"
)
ax1.legend()

#Rohdaten plotten (mit noise)
fig2, ax2 = plt.subplots(1, 1, layout="constrained")
ax2.plot(phi_m, A_m, "x", label="Amplituden mit Rauschen")
ax2.set(
    xlabel = r"$\phi / rad$",
    ylabel = r"$A / V$"
   # yscale = "linear"
)
ax2.legend()

#Rohdaten plotten (Photodetektor)
fig3, ax3 = plt.subplots(1, 1, layout="constrained")
ax3.plot(x_cm, U, "x", label="Photodetektor")
ax3.set(
    xlabel = r"$x / m$",
    ylabel = r"$U / V$",
    xlim = (0, 90),
    ylim = (0, 8)
)
ax3.legend()

plt.show()