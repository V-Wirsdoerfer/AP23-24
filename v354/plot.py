import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit 

#Daten generieren
t_p, U_p = np.genfromtxt("content/positive_Amplitude.txt", unpack=True)
t_n, U_n = np.genfromtxt("content/negative_Amplitude.txt", unpack=True)
f, U_C = np.genfromtxt("content/resonanzfrequenz.txt", unpack=True)

U_Res = U_C / 2.5
t1 = t_p * 1e6
t2 = t_n * 1e6
t_Ges = np.concatenate((np.array(t1), np.array(t2)))
U_Ges= np.concatenate((np.array((U_p)), np.array(U_n)))

#Definition für die curve fits
def exp1(t, a, mu, nu, eta):
    return a * np.e ** (-2*np.pi*mu*t) * np.cos(2*np.pi*nu*t + eta)

#plot zu a)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(t1, U_p, "rx", label="positive Amplituden")
ax.plot(t2, U_n, "bx", label="negative Amplituden")

#curve fit a)
params_n, cov_n = curve_fit(
    exp1,
    t_Ges,
    U_Ges,
    p0 = (1, 600000000, 3760, 0)
)
print(*params_n)
t = np.linspace(0, 500)
ax.plot(t, exp1(t, *params_n), label="curvefit")
ax.set(
    xlabel = "Zeit in Mikrosekunden",
    ylabel = "Spannung in Volt",
)

ax.legend()
#plt.show()

#plot zu c)
fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(f, U_C, "rx", label="Resonanzkurve")
ax.set(
    xlabel = "Frequenz in Hertz",
    ylabel = "Kondensatorspannung in Volt",
)

ax.legend()
plt.show()

