import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties import ufloat
from scipy.stats import sem

######################
# Konstanten definieren
######################

mu_0 = 1.256637062 * 10**-6  # in N/A^2
g = 9.81  # in N/kg
R_Spule = 0.109  # in m
d_Spule = 0.138  # in m
N_Spule = 195  # Anzahl der Windungen
M_Kugel = 0.15  # in kg
R_Kugel = 0.025  # in m
M_Masse = 1.4e-3  # in kg
r_offset = R_Kugel + 0.0139 - 0.0338  # in m, Abstand Mittelpunkt zu Ende Stabführung

# J_K = (2 * 0.15 * (0.025 ** 2)) / (5)
J_K = 2 / 5 * M_Kugel * R_Kugel**2
# b1 = 195 * ((mu_0) * (0.109 ** 2)) / ((0.109 ** 2) + (0.138/2) ** 2) ** 1.5


#################
# Daten generieren
#################

r_Grav, I_Grav = np.genfromtxt("./content/gravitation.txt", unpack=True)  # in cm und A
r_Grav /= 100  # in m umrechnen von cm
r_Grav += r_offset  # in m Gesamtabstand zum Mittelpunkt
I_Schwing, T_Schwing = np.genfromtxt(
    "./content/schwingung.txt", unpack=True
)  # in A und s
I_praz, T1_praz, T2_praz, T3_praz = np.genfromtxt(
    "./content/prazession.txt", unpack=True
)  # in A, s, s, s


######################
# Funktionen definieren
######################


def get_B(I):
    return (N_Spule * mu_0 * I * R_Spule**2) / (
        (R_Spule**2 + (d_Spule / 2) ** 2) ** (3 / 2)
    )


# def get_mu_Gravitation(I,r):
#    return (M_Masse * r * g) / (get_B(I))
#
# def get_mu_Schwingungsdauer(I, T):
#    return (4 * np.pi **2 * J_K) / (T**2 * get_B(I))
#
# def get_mu_Prazession(I, T):
#    L_Kugel = J_K * 2 * np.pi / T                   #L_K berechnen, da nicht direkt gegeben
#    return (2 * np.pi * L_Kugel) / (T * get_B(I))


####################
# Konstanten ausgeben
print("mu_0: ", mu_0)
# print("b1:", b1)
print("B(1A): ", get_B(1))
print("J_K: ", J_K)
print("r offset: ", r_offset)


###############################################################
# Plotten und aus Ausgleichsgerade magnetisches Moment berechnen
###############################################################


# Gravitation
############

# Daten
fig1, ax1 = plt.subplots(layout="constrained")
#ax1.plot(get_B(I_Grav) / M_Masse * g, r_Grav, "x", label="Datenpunkte")
ax1.plot(get_B(I_Grav), r_Grav, "x", label="Datenpunkte")

# linregress
#params1, cov1 = np.polyfit(get_B(I_Grav) / M_Masse * g, r_Grav, deg=1, cov=True)
params1, cov1 = np.polyfit(get_B(I_Grav), r_Grav, deg=1, cov=True)

err1 = np.sqrt(np.diag(cov1))
x = np.linspace(0.0029, 0.0041)
ax1.plot(x, x * params1[0] + params1[1], label="Ausgleichsgerade")

# settings
ax1.set(
    xlabel=r"$B / T$",
    ylabel=r"$r$/m",
    xlim=(min(x), max(x)),
)
ax1.legend()
fig1.savefig("build/Gravitation.pdf")

# Ausgabe Ergebnis
m_Grav = ufloat(params1[0], err1[0])
B_Grav = ufloat(params1[1], err1[1])
mu_Grav = M_Masse * g * m_Grav
print("\nmu_Grav, also Steigung: ", mu_Grav)
print("B_Grav: ", B_Grav)
print("mg*steigung, also mu: ", M_Masse*g*mu_Grav)


# Schwingung
###########

# Daten
fig2, ax2 = plt.subplots(layout="constrained")
#ax2.plot(
#    (T_Schwing**2) / (4 * J_K * np.pi**2),
#    1 / get_B(I_Schwing),
#    "x",
#    label="Datenpunkte",
#)
ax2.plot(1 / get_B(I_Schwing), T_Schwing**2, "x", label="Datenpunkte")


# linregress
#params2, cov2 = np.polyfit( (T_Schwing**2) / (4 * J_K * np.pi**2), 1 / get_B(I_Schwing), deg=1, cov=True)
params2, cov2 = np.polyfit(1 / get_B(I_Schwing), (T_Schwing**2), deg=1, cov=True)

err2 = np.sqrt(np.diag(cov2))
x = np.linspace(290, 1201)
ax2.plot(x, x * params2[0] + params2[1], label="Ausgleichsgerade")

# settings
ax2.set(
    xlabel=r"$\frac{1}{B} / T^{-1}$",
    ylabel=r"$T^2 / s^2$",
    xlim=(min(x), max(x)),
)
ax2.legend()
fig2.savefig("build/Schwingung.pdf")

# Ausgabe Ergebnis
m_Schwing = ufloat(params2[0], err2[0])
B_Schwing = ufloat(params2[1], err2[1])
mu_Schwing = (J_K*4 * np.pi**2)/(m_Schwing)
print("\nmu_Schwing: ", mu_Schwing)
print("B_Schwing: ", B_Schwing)
print("Steigung: ", m_Schwing)


# Präzession
###########

# Mittelwert der t's berechnen
T_praz = np.zeros(np.size(T1_praz))
T_praz_err = np.zeros(np.size(T1_praz))

for i in range(np.size(T1_praz)):
    t_arr = [T1_praz[i], T2_praz[i], T3_praz[i]]
    T_praz[i] = np.mean(t_arr)
    T_praz_err[i] = sem(t_arr)

T_praz_total = unp.uarray(T_praz, T_praz_err)

#x = get_B(I_praz) / (2 * np.pi * J_K * 2 * np.pi / T_praz_total)
x = get_B(I_praz)
x2 = (get_B(I_praz) * T_praz_total) / (4 * np.pi **2 * J_K)
y = 1 / T_praz_total

#print("\nT_praz_total: ", T_praz_total)
#print("\ny: ", y)
#print("\nx: ", x)
#print("\nx2: ", x2)
#print(x==x2)

fig3, ax3 = plt.subplots(layout="constrained")
ax3.errorbar(
    unp.nominal_values(x),
    unp.nominal_values(y),
    xerr=unp.std_devs(x),
    yerr=unp.std_devs(y),
    fmt="x",
    label="Datenpunkte",
)
#ax3.plot(get_B(I_praz) / (2 * np.pi * J_K * 2 * np.pi / T_praz), 1 / T_praz, "x", label="Datenpunkte normal" )

params3, cov3 = np.polyfit(unp.nominal_values(x),unp.nominal_values(y),deg=1,cov=True)
err3 = np.sqrt(np.diag(cov3))

lin_x = np.linspace(0.0012,0.006)
ax3.plot(lin_x, lin_x * params3[0] + params3[1], label="Ausgleichsgerade")

ax3.set(
    xlabel=r"$B / T$",
    ylabel=r"$\frac{1}{T} / s^{-1}$",
    xlim=(min(lin_x),max(lin_x))
)
ax3.legend()
fig3.savefig("build/Präzession.pdf")

#Ausgabe Ergebnisse
m_Praz= ufloat(params3[0], err3[0])
B_Praz  = ufloat(params3[1], err3[1])
mu_Praz = m_Praz*2 * np.pi * J_K * 2 * np.pi * 5.5

print("\nmu_Praz: ", mu_Praz)
print("B_Praz: ", B_Praz)
print("m Steigung:", m_Praz)