from turtle import bk
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp


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



### Grundgrößen berechnen ###



### Daten plotten ###

#Kupfer
fig_k, ax_k = plt.subplots(label="Kupfer")
ax_k.plot(I_Kupfer, UH_Kupfer, "x", label="Kupfer")
ax_k.set(
    xlabel="I/A",
    ylabel=r"U_H / V"
)
params_k, cov = np.polyfit(I_Kupfer, UH_Kupfer, deg=1, cov=True)
x = np.linspace(1,4)
ax_k.plot(x, params_k[0]*x + params_k[1], label = "ax + b")
ax_k.legend()
fig_k.savefig("./build/Kupfer.pdf")

#Silber
fig_s, ax_s = plt.subplots(label="Silber", layout = "constrained")
ax_s.plot(I_Silber, UH_Silber, "x", label="Silber")
ax_s.set(
    xlabel="I/A",
    ylabel=r"U_H / V"
)
params_s, cov = np.polyfit(I_Silber, UH_Silber, deg=1, cov=True)
x = np.linspace(1,4)
ax_s.plot(x, params_s[0]*x + params_s[1], label = f" ax + b \n a = {(params_s[0]*1e6):.3f} e-6")
ax_s.legend()
fig_s.savefig("./build/Silber.pdf")


