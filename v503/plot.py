import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.stats import sem


### Konstanten und Grundgrößen ###

g = 9.81            #in m/s^2
rho_Oel = 886       #in kg/m^3
rho_L = 1.1839      #in kg/m^3 bei T=25°C

### viskosität Luft bestimmen ###
T_Luft, eta_Luft = np.genfromtxt("data/rekonstruktion_eta.txt", unpack=True)
params, cov = np.polyfit(T_Luft, eta_Luft, deg=1, cov=True)
eta_Steigung = ufloat(params[0], np.sqrt(np.diag(cov))[0])
eta_0 = ufloat(params[1],np.sqrt(np.diag(cov))[1])
#print(params[0], np.sqrt(np.diag(cov))[0], eta_Steigung)
#x=np.linspace(15, 32)
#fig, ax = plt.subplots()
#plt.grid(True)
#ax.plot(x, params[0]*x + params[1])
#ax.plot(T_Luft, eta_Luft, "x")
#ax.set(xlim=(15,32))
#fig.savefig("build/eta.pdf")



### Funktionen und Formeln ###

def eta_Luft(T):
    return eta_Steigung * T + eta_0

def Ladung(T, v_ab, v_auf):
    q = 3 * np.pi * eta_Luft(T) * (v_ab + v_auf) * np.sqrt( (9* eta_Luft(T) * (v_ab - v_auf) ) / (4 * g * (rho_Oel - rho_L)) )
    return q

def Radius(T, v_ab, v_auf):
    r = np.sqrt( (9 * eta_Luft(T) * (v_ab - v_auf)) / (2 * g * (rho_Oel- rho_L)) )
    return r


