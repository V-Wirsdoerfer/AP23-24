from turtle import bk
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.stats import sem
from uncertainties import nominal_value


### Maaße und Grundgrößen ####

e0 = -1.60217662e-19
m0 = 9.10938356e-31
j = 1                       #in A/mm^2
h = 6.626e-34               #in Js

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

#Messmethode konstanter Probenstrom
I_Kupfer = 7                #in A  
I_Silber = 10               #in A   
I_Zink   = 7                #in A  

#Messmethode konstanter Spulenstrom
B_const = 0.576             #in T


### Daten importieren ###
I_Spule_Kupfer, U_Kupfer, B_Kupfer, UH_Kupfer = np.genfromtxt("./content/kupfer.txt", unpack = True)  #in A, V, mT, mV 
I_Spule_Silber, U_Silber, B_Silber, UH_Silber = np.genfromtxt("./content/silber.txt", unpack = True)  #in A, V, mT, mV
I_Spule_Zink, U_Zink, B_Zink, UH_Zink = np.genfromtxt("./content/zink.txt", unpack = True)            #in A, V, mT, mV
I_constB_auf, UH_constB_auf, I_constB_ent, UH_constB_ent = np.genfromtxt("./content/constB.txt", unpack = True) #aufladen und entladen


###in SI umrechnen
j *= 1e+6                      #in A/m^2

B_Kupfer *= 1e-3
UH_Kupfer *= 1e-3

B_Silber *= 1e-3
UH_Silber *= 1e-3

B_Zink *= 1e-3
UH_Zink *= 1e-3

UH_constB_auf *= 1e-3
UH_constB_ent *= 1e-3


#für Umpolung mit Vorzeichen berücksichtigen
I_constB_ent = -I_constB_ent
UH_constB_ent = -UH_constB_ent



### Funktionen ###
def Steigung(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[0], err[0])

def Ladungstraegerdichte(I, m, d):                 #Gleichung 6
    n = - (I) / (e0 * m * d)
    return n

def Ladungsträgerdichte_constB(m,d):
    n = - B_const / (e0 * m * d)
    return n

def mittlere_Flugzeit(L, n, Q, R):                  #Gleichung 5
    tau = (2 * m0 * L) / (e0**2 * n * Q * R)
    return tau

def mittlere_Driftgeschwindigkeit(UH, d, tau):      #Gleichung 1 und 2
    v_d = - (0.5 * e0 * UH * tau) / (m0 * d)
    #print("v_d: ", v_d)
    #print("np.sum/len:", np.sum(v_d)/len(v_d))
    #print("np.mean: ", np.mean(v_d))
    return np.mean(v_d)

def mittlere_Driftgeschwindigkeit_j(n):
    v_d = - j / (n * e0)
    return v_d

def Beweglichkeit(tau):                             #Gleichung 8
    µ = - (e0*tau)/(2*m0)
    return µ

def Fermi_Energie(n):
    E_f = h**2 / (2*m0) * pow(( ( (3*n) / (8*np.pi) )**2 ), (1/3))
    return E_f

def totalgeschwindigkeit(n):
    v = unp.sqrt((2*Fermi_Energie(n)) / (m0))
    return v

def mittlere_freie_Weglaenge(tau, v):               #Gleichung 7
    l = tau * abs(v)
    return l


### Daten plotten ###
#
#Orientierung an Gleichung 6
#Kupfer
fig_k, ax_k = plt.subplots(label="Kupfer", layout = "constrained")
plt.grid(True)
ax_k.plot(B_Kupfer, UH_Kupfer, "x", label="Kupfer")
ax_k.set(
    xlabel="B/T",
    ylabel=r"U_H / V",
    #xlim=[0,5]
)
params_k, cov_k = np.polyfit(B_Kupfer, UH_Kupfer, deg=1, cov=True)
m_k = Steigung(params_k,cov_k)
x = np.linspace(0,1.3)
ax_k.plot(x, params_k[0]*x + params_k[1], label = f"ax + b \na = {m_k}")
ax_k.legend()
fig_k.savefig("./build/Kupfer.pdf")


#Silber
fig_s, ax_s = plt.subplots(label="Silber", layout = "constrained")
plt.grid(True)
ax_s.plot(B_Silber, UH_Silber, "x", label="Silber")
ax_s.set(
    xlabel="B/T",
    ylabel=r"U_H / V",
    #xlim=[0,5]
)
params_s, cov_s = np.polyfit(B_Silber, UH_Silber, deg=1, cov=True)
m_s = Steigung(params_s,cov_s)
x = np.linspace(0,1.3)
ax_s.plot(x, params_s[0]*x + params_s[1], label = f" ax + b \n a = {m_s}")
ax_s.legend()
fig_s.savefig("./build/Silber.pdf")

#Zink
fig_z, ax_z = plt.subplots(label="Zink", layout = "constrained")
plt.grid(True)
ax_z.plot(B_Zink, UH_Zink, "x", label="Zink")
ax_z.set(
    xlabel="B/T",
    ylabel=r"U_H / V",
    #xlim=[0,4.8]
)
params_z, cov_z = np.polyfit(B_Zink, UH_Zink, deg=1, cov=True)
m_z = Steigung(params_z,cov_z)
x = np.linspace(0,1.1)
ax_z.plot(x, params_z[0]*x + params_z[1], label = f"ax + b \na = {m_z}")
ax_z.legend()
fig_z.savefig("./build/Zink.pdf")


### konstantes B-Feld ###
fig_constB, ax_constB = plt.subplots(layout="constrained")
plt.grid(True)
ax_constB.plot(I_constB_auf, UH_constB_auf, "rx", label="positiver Probenstrom")
ax_constB.plot(I_constB_ent, UH_constB_ent, "bx", label="negativer Probenstrom")

I_const = np.concatenate((I_constB_auf, I_constB_ent))
UH_const = np.concatenate((UH_constB_auf, UH_constB_ent))

params_constB, cov_constB = np.polyfit(I_const, UH_const, deg=1, cov=True)
m_constB = Steigung(params_constB, cov_constB)
x=np.linspace(-10.5,10.5)
ax_constB.plot(x, params_constB[0]*x + params_constB[1], "g", label=f" ax + b \n a = {m_constB}")

ax_constB.set(
    xlabel="I/A",
    ylabel="U_H/V",
)
ax_constB.legend()
fig_constB.savefig("./build/constB.pdf")

#
#
#
### Grundgrößen berechnen ###
#
#
#

### Ladungsträgerdichte berechnen mit Gleichung 6 ###

#Kupfer
n_Kupfer = Ladungstraegerdichte(I_Kupfer, m_k, Kupfer_Folie_d)
print("Ladungsträgerdichte n von Kupfer ist: ", n_Kupfer, "1/m^3")


#Silber
n_Silber_constI = Ladungstraegerdichte(I_Silber, m_s, Silber_Folie_d)
n_Silber_constB = Ladungsträgerdichte_constB(m_s, Silber_Folie_d)
print("Ladungsträgerdichte n von Silber bei konstanten Probenstrom ist: ", n_Silber_constI, "1/m^3")
print("Ladungsträgerdichte n von Silber bei konstanten B-Feld ist: ", n_Silber_constB, "1/m^3")
### hier auswählen mit welcher Methode weiter gerechnet werden soll ###
#n_Silber = n_Silber_constB
n_Silber = n_Silber_constI


#Zink
n_Zink = Ladungstraegerdichte(I_Zink, m_z, Zink_Folie_d)
print("Ladungsträgerdichte n von Zink ist: ", n_Zink, "1/m^3")


### Mittlereflugzeit tau berechnen mit Gleichung 5 ###

#Kupfer
tau_Kupfer = mittlere_Flugzeit(Kupfer_Spule_l, n_Kupfer, np.pi * (0.5*Kupfer_Draht_d)**2, Kupfer_Spule_R)
print("mittlere Flugzeit eines elektrons im Kupfer ist: ", tau_Kupfer, "s")

#Silber
tau_Silber = mittlere_Flugzeit(Silber_Spule_l, n_Silber, np.pi * (0.5*Silber_Draht_d)**2, Silber_Spule_R)
print("mittlere Flugzeit eines elektrons im Silber ist: ", tau_Silber, "s")

#Zink #### zu wenig infos
#tau_Zink = mittlere_Flugzeit(Zink_Spule_l, n_Zink, np.pi * (0.5*Zink_Draht_d)**2, Zink_Spule_R)
#print("mittlere Flugzeit eines elektrons im Zink ist: ", tau_Zink, "s")

#Tantal ### zu wenig infos
#tau_Tantal = mittlere_Flugzeit(Tantal_Spule_l, n_Tantal, np.pi * (0.5*Tantal_Draht_d)**2, Tantal_Spule_R)
#print("mittlere Flugzeit eines elektrons im Tantal ist: ", tau_Tantal, "s")
###


### mittlere Driftgeschwindigkeit berechnen mit Gleichung 1 und 2 ###

#Kupfer
#v_d_Kupfer = mittlere_Driftgeschwindigkeit(UH_Kupfer, Kupfer_Folie_d, tau_Kupfer)
#print("mittlere Driftgeschwindigkeit Kupfer: ", v_d_Kupfer, " m/s")

#Silber
#v_d_Silber = mittlere_Driftgeschwindigkeit(UH_Silber, Silber_Folie_d, tau_Silber)
#print("mittlere Driftgeschwindigkeit Silber: ", v_d_Silber, " m/s")

#Zink ## geht nicht, zu wenig Infos ##
#v_d_Zink = mittlere_Driftgeschwindigkeit(UH_Zink, Zink_Folie_d, tau_Zink)
#print("mittlere Driftgeschwindigkeit Zink: ", v_d_Zink, " m/s")

#Tantal ## geht nicht, zu wenig Infos ##
#v_d_Tantal = mittlere_Driftgeschwindigkeit(UH_Tantal, Tantal_Folie_d, tau_Tantal)
#print("mittlere Driftgeschwindigkeit Tantal: ", v_d_Tantal, " m/s")


### mittlere Driftgeschwindigkeit berechnen mit Gleichung 3 ###

#Kupfer
v_d_Kupfer = mittlere_Driftgeschwindigkeit_j(n_Kupfer)
print("mittlere Driftgeschwindigkeit Kupfer: ", v_d_Kupfer, " m/s")

#Silber
v_d_Silber = mittlere_Driftgeschwindigkeit_j(n_Silber)
print("mittlere Driftgeschwindigkeit Silber: ", v_d_Silber, " m/s")

#Zink
v_d_Zink = mittlere_Driftgeschwindigkeit_j(n_Zink)
print("mittlere Driftgeschwindigkeit Zink: ", v_d_Zink, " m/s")

#Tantal ## geht nicht, zu wenig Infos ##
#v_d_Tantal = mittlere_Driftgeschwindigkeit(UH_Tantal, Tantal_Folie_d, tau_Tantal)
#print("mittlere Driftgeschwindigkeit Tantal: ", v_d_Tantal, " m/s")


### Beqeglichkeit µ berechnen mit Gleichung 8 ###

#Kupfer
µ_Kupfer = Beweglichkeit(tau_Kupfer)
print("Beweglichkeit µ von Kupfer ist: ", µ_Kupfer)

#Silber
µ_Silber = Beweglichkeit(tau_Silber)
print("Beweglichkeit µ von Silber ist: ", µ_Silber)

#Zink ### nicht genug infos
#µ_Zink = Beweglichkeit(tau_Zink)
#print("Beweglichkeit µ von Zink ist: ", µ_Zink)


### totale Geschwindigkeit mit Abschätzung berechnen ###

#Kupfer
v_tot_Kupfer = totalgeschwindigkeit(n_Kupfer)
print("Totalgeschwindigkeit von Kupfer ist: ", v_tot_Kupfer)

#Silber
v_tot_Silber = totalgeschwindigkeit(n_Silber)
print("Totalgeschwindigkeit von Silber ist: ", v_tot_Silber)

#Zink
v_tot_Zink = totalgeschwindigkeit(n_Zink)
print("Totalgeschwindigkeit von Zink ist: ", v_tot_Zink)


### mittlere freie Weglänge l berechnen mit GLeichung 7 ###

#Kupfer
l_Kupfer = mittlere_freie_Weglaenge(tau_Kupfer, v_tot_Kupfer)
print("mittlere freie Weglänge l von Kupfer: ", l_Kupfer)

#Silber
l_Silber = mittlere_freie_Weglaenge(tau_Silber, v_tot_Silber)
print("mittlere freie Weglänge l von Silber: ", l_Silber)

#Zink ### nicht genügend infos 
#l_Zink = mittlere_freie_Weglaenge(tau_Zink, v_tot_Zink)
#print("mittlere freie Weglänge l: ", l_Zink)

