import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

c = const.c
e_0 = const.e
pi= const.pi
epsilon_0 = const.epsilon_0
n = 3.365
B_max = 429e-3  #B-Feld in Tesla

N_1 = 1.2e18
N_2 = 2.8e18 

d_0 = 5.11e-3   #Dicke in m
d_1 = 1.36e-3   #Dicke in m
d_2 = 1.296e-3  #Dicke in m

#Daten generieren

x, B = np.genfromtxt("content/Magnetfeld.txt", unpack=True)
lmbd, deg_0_rr, min_0_rr, deg_0_rb, min_0_rb = np.genfromtxt("content/hochrein.txt", unpack=True)    #reines GaAs ohne n-Dotierung
lmbd, deg_1_rr, min_1_rr, deg_1_rb, min_1_rb = np.genfromtxt("content/1_2e18.txt", unpack=True)      #leicht n-dotiertes GaAs
lmbd, deg_2_rr, min_2_rr, deg_2_rb, min_2_rb = np.genfromtxt("content/2_8e18.txt", unpack=True)      #stark n-dotiertes GaAs
lmbd *= 10e-7
B *= 10e-4

#Einehiten umrechnen

for i in [min_0_rb, min_0_rr, min_1_rb, min_1_rr, min_2_rb, min_2_rr]:
    i /= 60

#Genaue Winkel definieren (mit Minuten)

deg_0_exkt_rr = deg_0_rr + min_0_rr
deg_0_exkt_rb = deg_0_rb + min_0_rb
deg_1_exkt_rr = deg_1_rr + min_1_rr
deg_1_exkt_rb = deg_1_rb + min_1_rb
deg_2_exkt_rr = deg_2_rr + min_2_rr
deg_2_exkt_rb = deg_2_rb + min_2_rb

#Winkeldurchschnitt bestimmen

deg_0 = 0.5 * (deg_0_exkt_rb - deg_0_exkt_rr)
deg_1 = 0.5 * (deg_1_exkt_rb - deg_1_exkt_rr)
deg_2 = 0.5 * (deg_2_exkt_rb - deg_2_exkt_rr)

#Winkel in Radiant umrechnen

rad_0 = np.pi / 180 * deg_0
rad_1 = np.pi / 180 * deg_1
rad_2 = np.pi / 180 * deg_2

#Magnetfeldstärke plotten

fig1, ax1 = plt.subplots(layout="constrained")
ax1.errorbar(x, B, 0.5e-3, 0.5, fmt='rx', label="Messwerte")
ax1.legend()
ax1.set(
    xlabel=r"$z$ in mm",
    ylabel=r"$B$ in T"
)
ax1.grid()
fig1.savefig("Magnetfeld.pdf")

#Theta_Kr gegen lambda² auftragen

y_err = np.pi / (120 * 180)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.errorbar(lmbd**2, rad_0, xerr=None, yerr=y_err, fmt='rx', label="hochreines GaAs")
ax.errorbar(lmbd**2, rad_1, xerr=None, yerr=y_err, fmt='bx', label="1. dotiertes GaAs")
ax.errorbar(lmbd**2, rad_2, xerr=None, yerr=y_err, fmt='gx', label="2. dotiertes GaAs")
ax.set(
    xlabel = r"$\lambda^{2}$ in m²",
    ylabel = r"$\theta_{\mathrm{Kr}}$ in rad"
)
ax.grid()
ax.legend()
fig.savefig("Kristallwinkel.pdf")

#maskierte arrays definieren

lmbd_masked = np.delete(lmbd, (1, 7))
rad_2_masked = np.delete(rad_2, (1, 7))
rad_0_masked = np.delete(rad_0, (1, 7))

#Wellenlänge squared gegen Winkel plotten

fig2, ax2 = plt.subplots(1, 1, layout="constrained")
theta_frei_leicht = rad_1 / d_1 - rad_0 / d_0
theta_frei_stark = rad_2_masked / d_2 - rad_0_masked / d_0

ax2.errorbar(lmbd **2, rad_1 / d_1 - rad_0 / d_0, xerr=None, yerr=y_err, fmt='bx', label="Leicht dotiert")
ax2.errorbar(lmbd_masked **2, rad_2_masked / d_2 - rad_0_masked / d_0, xerr=None, yerr=y_err, fmt='gx', label="Stark dotiert")

ax2.set(
    xlabel=r"$\lambda$ in m ",
    ylabel=r"$\theta_{\mathrm{frei}}$ in rad/m"
)

#Ausgleichsrechnung stark dotiert

params_stark, cov_stark = np.polyfit(lmbd_masked **2, theta_frei_stark, deg=1, cov=True)
def f(lmbd_masked):
    return params_stark[0] * lmbd_masked**2 + params_stark[1]
ax2.plot(lmbd_masked**2, f(lmbd_masked), 'g',  label="Ausgleichsgerade stark dotiert")

#Ausgleichsrechnung leicht dotiert

params_leicht, cov_leicht = np.polyfit(lmbd **2, theta_frei_leicht, deg=1, cov=True)
def f(lmbd):
    return params_leicht[0] * lmbd**2 + params_leicht[1]
ax2.plot(lmbd**2, f(lmbd), 'b', label="Ausgleichsgerade leicht dotiert")

ax2.grid()
ax2.legend()
fig2.savefig("LinRegress.pdf")

#Effektive Masse berechnen

m_eff_leicht = np.sqrt(N_1 * B_max * e_0**3 / (8 * pi**2 * epsilon_0 * c**3 * params_leicht[0] * n))
m_eff_stark = np.sqrt(N_2 * B_max * e_0**3 / (8 * pi**2 * epsilon_0 * c**3 * params_stark[0] * n))

print("Maximale Magnetfeldstärke")
print("Effektive Masse über leicht dotierte Probe: ", m_eff_leicht)
print("Effektive Masse über stark dotierte Probe: ", m_eff_stark)





