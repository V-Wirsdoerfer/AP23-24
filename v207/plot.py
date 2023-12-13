from statistics import covariance
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
from scipy.stats import sem
#from scipy.optimize import curve_fitpyplot as plt

#Daten generieren
tg_ou1, tg_ou2, tg_uo1, tg_uo2 = np.genfromtxt("./content/gross_20C.txt", unpack="True")
tk_ou, tk_uo = np.genfromtxt("./content/klein_20C.txt", unpack="True")
T , td_ou1, td_ou2, td_uo1, td_uo2 = np.genfromtxt("./content/dynamisch.txt", unpack="True")

#Statisch klein plotten
K_kl = 0.07640
rho_kl = (4.9528) / (4/3 * np.pi * (1.576/2)**3)
rho_Fl = 0.99821
def s_kl(t):
    return K_kl * (rho_kl - rho_Fl) * t

t = np.linspace(0, 20)
fig, ax = plt.subplots(layout="constrained")
ax.plot = (t, s_kl(t), "rx", label="yeet")
ax.set(
    xlabel="Fallzeiten",
    ylabel="Viskosität", 
)
ax.legend()
plt.show()
