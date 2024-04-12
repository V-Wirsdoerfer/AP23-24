from turtle import bk
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.stats import sem


### Maaße und Grundgrößen ####

e0 = 1.60217662e-19
m0 = 9.10938356e-31

Tantal_Spule_l = 1.73       #in m 
Kupfer_Spule_l = 1.37       #in m
Silber_Spule_l = 1.73       #in m

Tantal_Spule_R = 4.865      #in Ohm
Kupfer_Spule_R = 2.745      #in Ohm
Silber_Spule_R = 0.642      #in Ohm

Tantal_Draht_d = ufloat(0.248e-3,   0.0005e-3)          #in m
Kupfer_Draht_d = ufloat(0.1e-3,     0.0005e-3)          #in m
Silber_Draht_d = ufloat(0.257e-3,   0.0005e-3)          #in m

Tantal_Folie_d = ufloat(0.037e-3,   0.005e-3)           #in m
Kupfer_Folie_d = ufloat(0.030e-3,   0.005e-3)           #in m
Silber_Folie_d = ufloat(0.030e-3,   0.005e-3)           #in m
Zink_Folie_d   = ufloat(0.0362e-3,  0.005e-3)           #in m



### Daten importieren ###
I_Kupfer, U_Kupfer, B_Kupfer, UH_Kupfer = np.genfromtxt("./content/kupfer.txt", unpack = True)  #in A, V, mT, mV 
I_Silber, U_Silber, B_Silber, UH_Silber = np.genfromtxt("./content/silber.txt", unpack = True)  #in A, V, mT, mV
I_Zink, U_Zink, B_Zink, UH_Zink = np.genfromtxt("./content/zink.txt", unpack = True)            #in A, V, mT, mV
I_constB_auf, UH_constB_auf, I_constB_ent, UH_constB_ent = np.genfromtxt("./content/constB.txt", unpack = True) #aufladen und entladen


#in SI umrechnen
B_Kupfer = B_Kupfer * 1e-3
UH_Kupfer = UH_Kupfer * 1e-3

B_Silber = B_Silber * 1e-3
UH_Silber = UH_Silber * 1e-3

B_Zink = B_Zink * 1e-3
UH_Zink = UH_Zink * 1e-3


#für Umpolung
I_constB_ent = -I_constB_ent
UH_constB_ent = -UH_constB_ent



### Funktionen ###
def Steigung(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[0], err[0])

def Ladungstraegerdichte (m, d):                #Gleichung 6
    n = - (1) / (e0 * m * d)
    return n

def mittlere_Flugzeit(L, n, Q, R):              #Gleichung 5
    tau = (2 * m0 * L) / (e0**2 * n * Q * R)
    return tau

def mittlere_Driftgeschwindigkeit(UH, d):       #Gleichung 1 und 2
    v_d = - (0.5 * e0 * UH ) / (m0 * d)
    print("np.mean:", np.mean(v_d))
    print("np.sem: ", sem(v_d))
    #v_d = ufloat(np.mean(v_d), np.std(v_d))
    return v_d


### Daten plotten ###
#
#Orientierung an Gleichung 6
#Kupfer
fig_k, ax_k = plt.subplots(label="Kupfer", layout = "constrained")
ax_k.plot(I_Kupfer * B_Kupfer, UH_Kupfer, "x", label="Kupfer")
ax_k.set(
    xlabel="IB/AT",
    ylabel=r"U_H / V",
    xlim=[0,5]
)
params_k, cov_k = np.polyfit(I_Kupfer * B_Kupfer, UH_Kupfer, deg=1, cov=True)
m_k = Steigung(params_k,cov_k)
x = np.linspace(0,5)
ax_k.plot(x, params_k[0]*x + params_k[1], label = f"ax + b \na = {m_k}")
ax_k.legend()
fig_k.savefig("./build/Kupfer.pdf")


#Silber
fig_s, ax_s = plt.subplots(label="Silber", layout = "constrained")
ax_s.plot(I_Silber * B_Silber, UH_Silber, "x", label="Silber")
ax_s.set(
    xlabel="IB/AT",
    ylabel=r"U_H / V",
    xlim=[0,5]
)
params_s, cov_s = np.polyfit(I_Silber * B_Silber, UH_Silber, deg=1, cov=True)
m_s = Steigung(params_s,cov_s)
x = np.linspace(0,5)
ax_s.plot(x, params_s[0]*x + params_s[1], label = f" ax + b \n a = {m_s}")
ax_s.legend()
fig_s.savefig("./build/Silber.pdf")

#Zink
fig_z, ax_z = plt.subplots(label="Zink", layout = "constrained")
ax_z.plot(I_Zink * B_Zink, UH_Zink, "x", label="Zink")
ax_z.set(
    xlabel="IB/AT",
    ylabel=r"U_H / V",
    xlim=[0,4.8]
)
params_z, cov_z = np.polyfit(I_Zink * B_Zink, UH_Zink, deg=1, cov=True)
m_z = Steigung(params_z,cov_z)
x = np.linspace(0,4.8)
ax_z.plot(x, params_z[0]*x + params_z[1], label = f"ax + b \na = {m_z}")
ax_z.legend()
fig_z.savefig("./build/Zink.pdf")



### Grundgrößen berechnen ###
#
#
### Ladungsträgerdichte berechnen mit Gleichung 6 ###

#Kupfer
n_Kupfer = Ladungstraegerdichte(m_k, Kupfer_Folie_d)
print("Die Ladungsträgerdichte n von Kupfer ist: ", n_Kupfer, "1/m^3")

#Silber
n_Silber = Ladungstraegerdichte(m_s, Silber_Folie_d)
print("Die Ladungsträgerdichte n von Silber ist: ", n_Silber, "1/m^3")

#Zink
n_Zink = Ladungstraegerdichte(m_z, Zink_Folie_d)
print("Die Ladungsträgerdichte n von Zink ist: ", n_Zink, "1/m^3")


### Mittlereflugzeit tau berechnen mit Gleichung 5 ###

#Kupfer
tau_Kupfer = mittlere_Flugzeit(Kupfer_Spule_l, n_Kupfer, np.pi * (0.5*Kupfer_Draht_d)**2, Kupfer_Spule_R)
print("die mittlere Flugzeit eines elektrons im Kupfer ist: ", tau_Kupfer, "s")

#Silber
tau_Silber = mittlere_Flugzeit(Silber_Spule_l, n_Silber, np.pi * (0.5*Silber_Draht_d)**2, Silber_Spule_R)
print("die mittlere Flugzeit eines elektrons im Silber ist: ", tau_Silber, "s")

### zu wenig infos
#Tantal
#tau_Tantal = mittlere_Flugzeit(Tantal_Spule_l, n_Tantal, np.pi * (0.5*Tantal_Draht_d)**2, Tantal_Spule_R)
#print("die mittlere Flugzeit eines elektrons im Tantal ist: ", tau_Tantal, "s")
###


### mittlere Driftgeschwindigkeit berechnen mit Gleichung 1 und 2 ###

#Kupfer
v_d_Kupfer = mittlere_Driftgeschwindigkeit(UH_Kupfer, Kupfer_Folie_d)
print("mittlere Driftgeschwindigkeit Kupfer: ", v_d_Kupfer, " m/s")
