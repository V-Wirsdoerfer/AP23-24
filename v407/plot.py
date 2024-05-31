import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat


#Daten generieren
Winkel, Intensität_senkrecht, drittel_senkrecht, Intensität_parallel, drittel_parallel, Fehler_senkrecht, Fehler_parallel = np.genfromtxt("Messwerte.txt", unpack=True)  #in °, A, A, A, A

I0 = ufloat(0.48e-3, 0.01e-3)   #in A


#Manche Werte durch 3 teilen wegen Faktor 3 auf Messgerät
drittel_senkrecht = np.asarray(drittel_senkrecht, dtype=bool)
drittel_parallel  = np.asarray(drittel_parallel,  dtype=bool)

Intensität_senkrecht[drittel_senkrecht] /= 3
Intensität_parallel[drittel_parallel]   /= 3

Fehler_senkrecht[drittel_senkrecht] /= 3
Fehler_parallel[drittel_parallel]   /= 3

#Größenordnung retten
Fehler_senkrecht  /= 100
Fehler_parallel   /= 100


#print("Intensität parallel:\n", np.round(Intensität_parallel, 7))
#print("Fehler parallel:\n", np.round(Fehler_parallel, 7))

fig1, ax1 = plt.subplots(layout="constrained")

ax1.errorbar(Winkel,  Intensität_senkrecht, xerr=0.5, yerr=Fehler_senkrecht, fmt="x", label="senkrecht")
ax1.errorbar(Winkel,  Intensität_parallel,  xerr=0.5, yerr=Fehler_parallel, fmt="x", label="parallel")
ax1.set(
    xlabel="Winkel / °",
    ylabel="Photostrom / A",
)
ax1.legend()

fig1.savefig("build/Messdaten.pdf")

    
#wurzel I/I0
fig2, ax2 = plt.subplots(layout="constrained")

ax2.errorbar(Winkel,  np.sqrt(Intensität_senkrecht/unp.nominal_values(I0)),xerr=0.5, yerr=np.sqrt(Fehler_senkrecht/unp.nominal_values(I0)) , fmt="x", label="senkrecht")
ax2.errorbar(Winkel,  np.sqrt(Intensität_parallel/unp.nominal_values(I0)), xerr=0.5, yerr=np.sqrt(Fehler_parallel/unp.nominal_values(I0))  , fmt="x", label="parallel")
ax2.set(
    xlabel="Winkel / °",
    ylabel=r"$\sqrt{\frac{I}{I_0}}$",
)
ax2.legend(loc="upper left")

fig2.savefig("build/I0.pdf")

#Brewsterwinkel bestimmen
Brewster_Intensität = min(Intensität_parallel[0:-3])
Brewster_Winkel = Winkel[Intensität_parallel == Brewster_Intensität]

print("Intensität beim Brewsterwinkel: ", Brewster_Intensität)
print("Brewsterwinkel: ", Brewster_Winkel)


### Brechungsindex berechnen
E_r_senkrecht = unp.uarray(np.sqrt(Intensität_senkrecht), np.sqrt(Fehler_senkrecht))
E_r_parallel = unp.uarray(np.sqrt(Intensität_parallel), np.sqrt(Fehler_parallel))
E_einfall = unp.sqrt(I0)


def Brechungsindex_senkrecht(): 
    alpha = 2*np.pi * Winkel/360    #alpha in Radiant 
    n = unp.sqrt( (- 2 * E_einfall * E_r_senkrecht * np.cos(2*alpha) + E_r_senkrecht**2 + E_einfall**2) / (2 * E_einfall * E_r_senkrecht + E_r_senkrecht**2 + E_einfall**2) )
    return n

def Brechungsindex_parallel(): 
    alpha = 2*np.pi * Winkel/360    #alpha in Radiant
    n = 1/(2 * np.cos(alpha)**2 * (E_r_parallel - E_einfall)**2 ) + unp.sqrt(1 / (4 * np.cos(alpha)**4 * (E_r_parallel - E_einfall)**4) - np.tan(alpha)**2 * ( (E_r_parallel + E_einfall) / (E_r_parallel - E_einfall) )**2 )
    return n 

print("Brechungsindex: ", len(Brechungsindex_parallel()) )
print("len Winkel: ", len(Winkel))


### Brechungsindex einsehen
fig3, ax3 = plt.subplots(2,2,layout="constrained")

ax3 = np.ravel(ax3)
ax3[0].errorbar(Winkel, unp.nominal_values(Brechungsindex_parallel()), yerr=unp.std_devs(Brechungsindex_parallel()), fmt="x", label="Brechungsindex \nparallel")
ax3[1].errorbar(Winkel, np.log(unp.nominal_values(Brechungsindex_parallel())), yerr=np.log(unp.std_devs(Brechungsindex_parallel())), fmt="x", label="Brechungsindex \nparallel \nlogarithmiert")
ax3[2].errorbar(Winkel, unp.nominal_values(Brechungsindex_senkrecht()), yerr=unp.std_devs(Brechungsindex_senkrecht()), fmt="x", label="Brechungsindex \nsenkrecht")
ax3[3].errorbar(Winkel, np.log(unp.nominal_values(Brechungsindex_senkrecht())), yerr=abs(np.log(unp.std_devs(Brechungsindex_senkrecht()))), fmt="x", label="Brechungsindex \nsenkrecht \nlogarithmiert")


ax3[0].set(
    xlabel = "Winkel/°",
    ylabel = "n",
    yscale = "log"
)

ax3[2].set(
    xlabel = "Winkel/°",
    ylabel = "n",
    yscale = "log"
)

for i in range(len(ax3)):
    print(i)
    ax3[i].legend(loc = "upper left")

fig3.savefig("build/Brechungsindex.pdf")
