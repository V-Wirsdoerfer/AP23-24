import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from uncertainties import ufloat
from scipy.optimize import root_scalar
import uncertainties.unumpy as unp

#Daten generieren
t_p, U_p = np.genfromtxt("content/positive_Amplitude.txt", unpack=True)
t_n, U_n = np.genfromtxt("content/negative_Amplitude.txt", unpack=True)
f, U_C = np.genfromtxt("content/resonanzfrequenz.txt", unpack=True)

U_Res = U_C / 2.5
t1 = t_p * 1e6
t2 = t_n * 1e6
t_Ges = np.concatenate((np.array(t_p), np.array(t_n)))
U_Ges= np.concatenate((np.array((U_p)), np.array(U_n)))

#Grundsätzliche Größen
L = ufloat(0.01011, 0.00003)
C = ufloat(2.093e-9, 0.003e-9)
R_1 = ufloat(48.1,0.1)

#Definition für die curve fits
def exp1(t, a, mu, nu, eta):                                                            #für a)
    return a * np.e ** (-2*np.pi*mu*t) * np.cos(2*np.pi*nu*t + eta)
def exp2(x, a, sigma, mu):                                                              #für c)
    return a / (sigma * np.sqrt(2*np.pi)) * np.e ** (-0.5 * ((x-mu) / sigma)**2) + 1


#plot zu a)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(t_p, U_p, "rx", label="positive Amplituden")
ax.plot(t_n, U_n, "bx", label="negative Amplituden")

#curve fit a)
params_n, cov_n = curve_fit(
    exp1,
    t_Ges,
    U_Ges,
    p0 = (4, 600, 3760, 0)
)
print("parameter zu A:", *params_n)
t = np.linspace(0, 500 * 1e-6, 10000)
ax.plot(t, exp1(t, *params_n), label="curvefit")
ax.set(
    xlabel = "Zeit in Mikrosekunden",
    ylabel = "Spannung in Volt",
)

ax.legend()
fig.savefig("build/a.pdf")


#Resonanzfrequenz berechnen b)
R_ap = unp.sqrt(4*L/C)
print("Aperiodische Dämpfung bei Widerstand R_ap: ", R_ap)


#plot zu  c) (halblogarithmisch)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(f, U_Res, "rx", label="Resonanzkurve")
ax.set(
    xlabel = "Frequenz in Hertz",
    ylabel = "Kondensatorspannung in Volt",
    yscale = "log"
)

ax.legend()
fig.savefig("build/c_log.pdf")

#plot zu c) (linear)
params_c, cov_c = curve_fit(
    exp2,
    f,
    U_Res,
    p0 = (3.5, 5000, 31000)
)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(f, U_Res, "rx", label="Resonanzkurve")
x = np.linspace(0, 50000)
ax.plot(x, exp2(x, *params_c), label="curvefit")
ax.plot(x, 1 / np.sqrt(2) * max(U_Res) +x-x, "_", label="horizontaler Vincent")

#Schnittpunkt
w1 = root_scalar(lambda x: exp2(x, *params_c) - max(U_Res) / np.sqrt(2), bracket = [0, 35000], method = "bisect")
print("w1: ", w1.root)
w2 = root_scalar(lambda x: exp2(x, *params_c) - max(U_Res) / np.sqrt(2), bracket = [30000, 45000], method = "bisect")
print("w2: ", w2.root)
w0 = 1 / unp.sqrt(L*C)
q = w0 / (w2.root - w1.root)
print("q exp: ", q)
ax.set(
    xlabel = "Frequenz in Hertz",
    ylabel = "Kondensatorspannung in Volt",
)

err_c = np.sqrt(np.diag(cov_c))

print("params_c: ", params_c, "\nerr_c:", err_c)
ax.legend()
fig.savefig("build/c_lin.pdf")
#plt.show()

#Dämpfungswiderstand berechnen
error_params = np.sqrt(np.diag(cov_n))
a = ufloat(params_n[0], error_params[0])
mu = ufloat(params_n[1], error_params[1])
nu = ufloat(params_n[2], error_params[2])
eta = ufloat(params_n[3], error_params[3])
R = 4 * np.pi * L * mu
T = 2*L / R

w1_theo = R/(2*L) + unp.sqrt(R**2/(4*L**2) + 1/(L*C))
w2_theo = -R/(2*L) + unp.sqrt(R**2/(4*L**2) + 1/(L*C))
q_theo = 1/(w0 *R*C) 
print("w1 und w2 theo:", w1_theo, w2_theo)
print("q theoretisch berechnet:", q_theo)
print("Dämpfungswiderstand: ", R)
print("Abklingdauer: ", T)
#print("Parameter a: ", a)
#print("Parameter mu: ", mu)
#print("Parameter nu: ", nu)
#print("Parameter eta: ", eta)

#Güte q des Schwingkreises berechnen


print("R diff : ", R-R_1)