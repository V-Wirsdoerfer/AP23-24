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



### Daten generieren ###

auf1, ab1, T1 = np.genfromtxt("./data/drop1.txt", unpack=True)      #in s, s, °C
auf4, ab4, T4 = np.genfromtxt("./data/drop4.txt", unpack=True)      #in s, s, °C
auf6, ab6, T6 = np.genfromtxt("./data/drop6.txt", unpack=True)      #in s, s, °C
auf7, ab7, T7 = np.genfromtxt("./data/drop7.txt", unpack=True)      #in s, s, °C
auf8, ab8, T8 = np.genfromtxt("./data/drop8.txt", unpack=True)      #in s, s, °C
auf9, ab9, T9 = np.genfromtxt("./data/drop9.txt", unpack=True)      #in s, s, °C
auf10, ab10, T10 = np.genfromtxt("./data/drop10.txt", unpack=True)  #in s, s, °C    
auf11, ab11, T11 = np.genfromtxt("./data/drop11.txt", unpack=True)  #in s, s, °C    
auf12, ab12, T12 = np.genfromtxt("./data/drop12.txt", unpack=True)  #in s, s, °C    
auf13, ab13, T13 = np.genfromtxt("./data/drop13.txt", unpack=True)  #in s, s, °C    
auf14, ab14, T14 = np.genfromtxt("./data/drop14.txt", unpack=True)  #in s, s, °C    
auf15, ab15, T15 = np.genfromtxt("./data/drop15.txt", unpack=True)  #in s, s, °C    

### Mittelwerte berechnen ###

#auf = unp.array(12)
#print(auf)
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf4 = ufloat(np.mean(auf4), sem(auf4))
auf6 = ufloat(np.mean(auf6), sem(auf6))
auf7 = ufloat(np.mean(auf7), sem(auf7))
auf8 = ufloat(np.mean(auf8), sem(auf8))
auf9 = ufloat(np.mean(auf9), sem(auf9))
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf1 = ufloat(np.mean(auf1), sem(auf1))
auf1 = ufloat(np.mean(auf1), sem(auf1))




















### Funktionen und Formeln ###

def eta_Luft(T):
    return eta_Steigung * T + eta_0

def Ladung(T, v_ab, v_auf):
    q = 3 * np.pi * eta_Luft(T) * (v_ab + v_auf) * unp.sqrt( (9* eta_Luft(T) * (v_ab - v_auf) ) / (4 * g * (rho_Oel - rho_L)) )
    return q

def Radius(T, v_ab, v_auf):
    r = unp.sqrt( (9 * eta_Luft(T) * (v_ab - v_auf)) / (4 * g * (rho_Oel- rho_L)) )
    return r


