import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
from uncertainties import ufloat
from scipy.optimize import curve_fit


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
hilfs_Intensität_senkrecht = unp.uarray(Intensität_senkrecht, Fehler_senkrecht)
hilfs_Intensität_parallel = unp.uarray(Intensität_parallel, Fehler_parallel)

fig2, ax2 = plt.subplots(layout="constrained")

ax2.errorbar(Winkel,  unp.nominal_values(unp.sqrt(hilfs_Intensität_senkrecht/I0)),xerr=0.5, yerr=unp.std_devs(unp.sqrt(hilfs_Intensität_senkrecht/I0)) , fmt="x", label="senkrecht")
ax2.errorbar(Winkel,  unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)), xerr=0.5, yerr=unp.std_devs(unp.sqrt(hilfs_Intensität_parallel/I0))  , fmt="x", label="parallel")
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
    n = (1)/(2 * np.cos(alpha)**2 * (E_r_parallel - E_einfall)**2 ) + unp.sqrt(1 / (4 * np.cos(alpha)**4 * (E_r_parallel - E_einfall)**4) - np.tan(alpha)**2 * ( (E_r_parallel + E_einfall) / (E_r_parallel - E_einfall) )**2 )
    return n 

print("Anzahl Brechungsindizes: ", len(Brechungsindex_parallel()) )
print("Anzahl Winkel: ", len(Winkel))

n_p = Brechungsindex_parallel()
log_n_p = unp.log(n_p,10)
log_log_n_p = unp.log(log_n_p)
sqrt_log_n_p = unp.pow(log_n_p,1/2)


n_s = Brechungsindex_senkrecht()
log_n_s = unp.log(n_s)
#log_log_n_s = unp.log(log_n_s)     #log_n_s wird negativ
sqrt_n_s = unp.sqrt(n_s)
sqrt_log_n_s = unp.pow(abs(log_n_s),1/2)


### Brechungsindex einsehen parallel
fig3, ax3 = plt.subplot_mosaic([["a","b"], ["c","d"]], layout="constrained")

for loc, val, fmt in zip( ("acd"), (n_p, log_n_p, sqrt_log_n_p), ("xb", "xr", "xg") ):
    ax3[loc].errorbar(Winkel, unp.nominal_values(val), yerr=unp.std_devs(val), fmt=fmt, label="Brechungsindex \nparallel")
ax3["a"].set(
    xlabel = "Winkel/°",
    ylabel = r"$n_\parallel$",
    yscale = "log"
)
ax3["c"].set(
    xlabel = "Winkel/°",
    ylabel = r"$\ln{(n_\parallel)}$",
)
ax3["d"].set(
    xlabel = "Winkel/°",
    ylabel = r"$\sqrt{-\ln{(n_\parallel)}}$i",
)

ax3["b"].set_axis_off()
for label, fmt in zip(("\nBrechungsindex parallel \nunverändert","\nBrechungsindex parallel \nlogarithmiert ","\nBrechungsindex parallel \nlogarithmiert \nwurzel gezogen"), ("brg")):
    ax3["b"].plot(1,1, color=fmt, label=label)
ax3["b"].legend()

fig3.savefig("build/Brechungsindex_parallel.pdf")


### Brechungsindex orthogonal/senkrecht polarisiert
fig4, ax4 = plt.subplot_mosaic([["left", "legend"], ["lower left", "lower right"]], layout="constrained")

ax4["left"].errorbar(Winkel, unp.nominal_values(n_s), yerr=unp.std_devs(n_s), fmt="bx", label="Brechungsindex \nsenkrecht")
ax4["lower left"].errorbar(Winkel, unp.nominal_values(log_n_s), yerr=unp.std_devs(log_n_s), fmt="rx", label="Brechungsindex \nsenkrecht logarithmiert")
ax4["lower right"].errorbar(Winkel, unp.nominal_values(sqrt_log_n_s), yerr=abs(unp.std_devs(sqrt_log_n_s)), fmt="gx", label="Brechungsindex \nsenkrecht \nlogarithmiert gewurzelt")

ax4["legend"].plot(1,1,color ="b",label="\nBrechungsindex senkrecht")
ax4["legend"].plot(1,1,color ="r",label="\nBrechungsindex senkrecht \nlogarithmiert")
ax4["legend"].plot(1,1,color ="g",label="\nBrechungsindex senkrecht \nlogarithmiert \nwurzel gezogen")

ax4["legend"].set_axis_off()
ax4["legend"].legend()

ax4["left"].set(
    xlabel = "Winkel/°",
    ylabel = r"$n_\bot$",
    yscale = "log"
)
ax4["lower left"].set(
    xlabel = "Winkel/°",
    ylabel = r"$\ln{(n_\bot)}$",
)
ax4["lower right"].set(
    xlabel = "Winkel/°",
    ylabel = r"$\sqrt{-\ln{(n_\bot)}}$i",
)

fig4.savefig("build/Brechungsindex_orthogonal.pdf")



### Fit von Gleichung 18 und 21 aus Anleitung an Graphen

def E_r_orth(alpha, n, a):     #Gleichung 18 / sqrt(I0) also /E_einfall
    '''
    Gleichung 18
    Gibt die Amplitude der Reflektierten Strahlung zurück
    nimmt nur Werte ohne Fehler
    Gibt Werte mit Fehler zurück
    '''
    alpha = 2*np.pi * alpha / 360       #alpha in Radiant umrechnen
    E =  (np.sqrt(n**2 - (np.sin(alpha))**2 ) - np.cos(alpha) )**2 / (n**2 -1) + a
    return E

def E_r_parallel(alpha, n, a): #Gleichung 21 / sqrt(I0) also /E_einfall
    '''
    Gleichung 21
    Gibt die Amplitude der Reflektierten Strahlung zurück
    nimmt nur Werte ohne Fehler
    Gibt Werte mit Fehler zurück
    '''
    alpha = 2*np.pi * alpha / 360       #alpha in Radiant umrechnen
    E = (n**2 * np.cos(alpha) - np.sqrt(n**2 - (np.sin(alpha))**2) ) / ( n**2 * np.cos(alpha) + np.sqrt(n**2 - (np.sin(alpha))**2) ) + a
    return E

params_senkrecht, cov_senkrecht = curve_fit(E_r_orth, Winkel, unp.nominal_values(unp.sqrt(hilfs_Intensität_senkrecht/I0)), p0=[7,-0.37])

x_alpha = np.linspace(0, 90,10000)
ax2.plot(x_alpha, E_r_orth(x_alpha, 7, -0.37), label="geratener Plot")
ax2.plot(x_alpha, E_r_orth(x_alpha, *params_senkrecht), label ="fit")


params_parallel, cov_parallel = curve_fit(E_r_parallel, Winkel, unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)), p0=[7, -0.5])
params_parallel2, cov_parallel2 = curve_fit(E_r_parallel, Winkel[Winkel <=74], unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)[Winkel<=74]), p0=[7, -0.5])

x_alpha = np.linspace(0, 90,10000)
x_alpha2 = np.linspace(0, 74,10000)
ax2.plot(x_alpha, E_r_parallel(x_alpha, 12, -0.54), label="geratener Plot")
ax2.plot(x_alpha, E_r_parallel(x_alpha, *params_parallel), label ="fit")
ax2.plot(x_alpha2, E_r_parallel(x_alpha2, *params_parallel2), label ="fit")


ax2.legend()
ax2.set(ylim=[0,0.7])
fig2.savefig("build/I0_fit.pdf")

print("senkrecht, n, a, b:", *params_senkrecht)
print("parallel, n, a, b:", *params_parallel)

