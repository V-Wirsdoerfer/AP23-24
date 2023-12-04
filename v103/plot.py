from pyexpat.errors import XML_ERROR_FINISHED
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat

# Daten generieren
x_K_e, Dx_K_e = np.genfromtxt("./content/Kreis_einseitig.txt", unpack=True)
x_K_b_n, Dx_K_b_n = np.genfromtxt("./content/Kreis_beidseitig_nah.txt", unpack=True)
x_K_b_f, Dx_K_b_f = np.genfromtxt("./content/Kreis_beidseitig_fern.txt", unpack=True)
x_Q_e, Dx_Q_e = np.genfromtxt("./content/Quadrat_einseitig.txt", unpack=True)
x_Q_b_n, Dx_Q_b_n = np.genfromtxt("./content/Quadrat_beidseitig_nah.txt", unpack=True)
x_Q_b_f, Dx_Q_b_f = np.genfromtxt("./content/Quadrat_beidseitig_fern.txt", unpack=True)


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
def error_cov(cov):
    return np.sqrt(np.diag(cov))

def print_Steigung(params):
    print("Steigung F/(48 E I) der Ausgleichsgrade: ", params[0])
    
def get_Steigung(params, cov):
    return ufloat(params[0], error_cov(cov)[0])

def print_E_Modul_beidseitig(Name, F, I, Steigung):
    E = F/(48*I*Steigung)
    print("Der Elastizitätsmodul von ", Name, " beträgt: ", E)

def print_E_Modul_einseitig(Name, F, I, Steigung):
    E = F/(2*I*Steigung)
    print("Der Elastizitätsmodul von ", Name, " beträgt: ", E)


#def Dx_beidseitig_nah(x, E, L, F, I):
#    return (F / 48 * E * I) * (3 * (L**2) * x - 4 * x**3)
#
#
#def Dx_beidseitig_fern(x, E, L, F, I):
#    return (F / 48 * E * I) * (4 * x**3 - 12 * L * x**2 + 9 * (L**2) * x - L**3)
#
#
#def Dx_einseitig(x, E, L, F, I):
#    return (F / 2 * E * I) * (L * x**2 - (x**3) / 3)


# linearisieren
# einseitig linearisieren
lin_x_K_e = L_K_e.nominal_value * x_K_e**2 - (x_K_e**3) / 3
lin_x_Q_e = L_Q_e.nominal_value * x_Q_e**2 - (x_Q_e**3) / 3


# beidseitig linearisieren
lin_x_K_b_n = 3 * (L_K_b.nominal_value**2) * x_K_b_n - 4 * (x_K_b_n**3)
lin_x_K_b_f = (
    4 * (x_K_b_f**3)
    - 12 * L_K_b.nominal_value * (x_K_b_f**2)
    + 9 * (L_K_b.nominal_value**2) * x_K_b_f
    - L_K_b.nominal_value**3
)

lin_x_Q_b_n = 3 * (L_Q_b.nominal_value**2) * x_Q_b_n - 4 * (x_Q_b_n**3)
lin_x_Q_b_f = (
    4 * (x_Q_b_f**3)
    - 12 * L_Q_b.nominal_value * (x_Q_b_f**2)
    + 9 * (L_Q_b.nominal_value**2) * x_Q_b_f
    - L_Q_b.nominal_value**3
)


# ausgleichsgeraden erstellen
params_K_e, cov_K_e = np.polyfit(lin_x_K_e, Dx_K_e, deg=1, cov=True)
params_K_b_n, cov_K_b_n = np.polyfit(lin_x_K_b_n, Dx_K_b_n, deg=1, cov=True)
params_K_b_f, cov_K_b_f = np.polyfit(lin_x_K_b_f, Dx_K_b_f, deg=1, cov=True)

params_Q_e, cov_Q_e = np.polyfit(lin_x_Q_e, Dx_Q_e, deg=1, cov=True)
params_Q_b_n, cov_Q_b_n = np.polyfit(lin_x_Q_b_n, Dx_Q_b_n, deg=1, cov=True)
params_Q_b_f, cov_Q_b_f = np.polyfit(lin_x_Q_b_f, Dx_Q_b_f, deg=1, cov=True)


# plots erstellen
# kreisförmig einseitig
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(lin_x_K_e, Dx_K_e, "rx", label="kreisförmig, einseitig")
x = np.linspace(0, 0.06)
ax1.plot(x, params_K_e[0] * x + params_K_e[1], label="lineare Rergression")
ax1.set(
    xlabel=r"$Lx^2 - \frac{x^3}{3} /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.06),
    ylim=(0, 0.0025),
)

ax1.legend()
fig1.savefig("./build/kreis_e.pdf")

# kreisförmig beidseitig nah
fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(
    lin_x_K_b_n,
    Dx_K_b_n,
    "gx",
    label=r"kreisförmig, beidseitig: $0 < x < \frac{1}{2} L$",
)
x = np.linspace(0, 0.18)
ax2.plot(x, params_K_b_n[0] * x + params_K_b_n[1], label="lineare Regression")
ax2.set(
    xlabel=r"$3L^2x - 4x^3 /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.18),
    ylim=(params_K_b_n[1], 0.00025),
)
ax2.legend()
fig2.savefig("./build/kreis_b_n.pdf")

# kreisförmig beidseitig fern
fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(
    lin_x_K_b_f,
    Dx_K_b_f,
    "rx",
    label=r"kreisförmig, beidseitig: $\frac{1}{2} L < x < L$",
)
x = x  # x Werte von k_b_n übernommen
ax3.plot(x, params_K_b_f[0] * x + params_K_b_f[1], label="lineare Regression")
ax3.set(
    xlabel=r"$4 x ^3 -12 L x ^2 + 9L^2 x - L^3 /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.18),
    ylim=(0.000075, 0.000275),
)
ax3.legend()
fig3.savefig("./build/kreis_b_f.pdf")


# quadratisch einseitig
fig4, ax4 = plt.subplots(layout="constrained")
ax4.plot(lin_x_Q_e, Dx_Q_e, "rx", label="quadratisch, einseitig")
x = np.linspace(0, 0.07)
ax4.plot(x, params_Q_e[0] * x + params_Q_e[1], "g", label="lineare Regression")
ax4.set(
    xlabel=r"$3L^2x - 4x^3 /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.07),
    ylim=(0, 0.00175),
)
ax4.legend()
fig4.savefig("./build/quadrat_e.pdf")

# quadratisch beidseitig nah
fig5, ax5 = plt.subplots(layout="constrained")
ax5.plot(
    lin_x_Q_b_n,
    Dx_Q_b_n,
    "bx",
    label=r"quadratisch, beidseitig: $0 < x < \frac{1}{2} L$",
)
x = np.linspace(0, 0.18)
ax5.plot(x, params_Q_b_n[0] * x + params_Q_b_n[1], "r", label="lineare Regression")
ax5.set(
    xlabel=r"$3L^2x - 4x^3 /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.18),
    ylim=(params_Q_b_n[1], 0.00016),
)
ax5.legend()
fig5.savefig("./build/quadrat_b_n.pdf")

# quadratisch beidseitig fern
fig6, ax6 = plt.subplots(layout="constrained")
ax6.plot(
    lin_x_Q_b_f,
    Dx_Q_b_f,
    "rx",
    label=r"quadratisch, beidseitig $\frac{1}{2} L < x < L$",
)
x = x  # x von Q_b_n wird übernommen
ax6.plot(x, params_Q_b_f[0] * x + params_Q_b_f[1], "g", label="lineare Regression")
ax6.set(
    xlabel=r"$4 x ^3 -12 L x ^2 + 9L^2 x - L^3 /$ m³",
    ylabel=r"$D(x)/$ m",
    xlim=(0, 0.18),
    #ylim=(0, 0.00017),
)
ax6.legend()
fig6.savefig("./build/quadrat_b_f.pdf")


#print("Steigung von Kreis einseitig: ", get_Steigung(params_K_e, cov_K_e))
#print("Steigung von Kreis beidseitig nah: ", get_Steigung(params_K_b_n, cov_K_b_n))
#print("Steigung von Kreis beidseitig fern: ", get_Steigung(params_K_b_f, cov_K_b_f))
#print("Steigung von Quadrat einseitig: ", get_Steigung(params_Q_e, cov_Q_e))
#print("Steigung von Quadrat beidseitig nah: ", get_Steigung(params_Q_b_n, cov_Q_b_n))
#print("Steigung von Quadrat beidseitig fern: ", get_Steigung(params_Q_b_f, cov_Q_b_f))




#print_E_Modul_einseitig("Kreis einseitig", F_einseitig, I_K, get_Steigung(params_K_e, cov_K_e))
#print_E_Modul_beidseitig("Kreis beidseitig nah", F_beidseitig, I_K, get_Steigung(params_K_b_n, cov_K_b_n))
#print_E_Modul_beidseitig("Kreis beidseitig fern", F_beidseitig, I_K, get_Steigung(params_K_b_f, cov_K_b_f))
#print_E_Modul_einseitig("Quadrat einseitig", F_einseitig, I_Q, get_Steigung(params_Q_e, cov_Q_e))
#print_E_Modul_beidseitig("Quadrat beidseitig nah", F_beidseitig, I_Q, get_Steigung(params_Q_b_n, cov_Q_b_n))
#print_E_Modul_beidseitig("Quadrat beidseitig fern", F_beidseitig, I_Q, get_Steigung(params_Q_b_f, cov_Q_b_f))

#Dichte bestimmen
L_K = ufloat(0.59,0.005)
L_Q = ufloat(0.6,0.005)
D_K = ufloat(0.01,0.000005)
D_Q = ufloat(0.01,0.000005)
m_K = ufloat(0.4123,0.00005)
m_Q = ufloat(0.5361,0.00005)

rho_K = m_K / (L_K * (0.5*D_K)**2 * np.pi) 
rho_Q = m_Q / (L_Q * (D_Q)**2)

print("Dichte kreisförmig: ", rho_K, "\nDichte quadratisch: ", rho_Q)

#Werte ausgeben

