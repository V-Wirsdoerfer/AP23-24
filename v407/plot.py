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
Brewster_Winkel = ufloat(Winkel[Intensität_parallel == Brewster_Intensität], 0.5)

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
    E = a*abs((np.sqrt(n**2 - (np.sin(alpha))**2 ) - np.cos(alpha) )**2 / (n**2 -1)) 
    return E

def E_r_parallel(alpha, n, a): #Gleichung 21 / sqrt(I0) also /E_einfall
    '''
    Gleichung 21
    Gibt die Amplitude der Reflektierten Strahlung zurück
    nimmt nur Werte ohne Fehler
    Gibt Werte mit Fehler zurück
    '''
    alpha = 2*np.pi * alpha / 360       #alpha in Radiant umrechnen
    E = a*abs((n**2 * np.cos(alpha) - np.sqrt(n**2 - (np.sin(alpha))**2) ) / ( n**2 * np.cos(alpha) + np.sqrt(n**2 - (np.sin(alpha))**2) ))
    return E

params_senkrecht, cov_senkrecht = curve_fit(E_r_orth, Winkel, unp.nominal_values(unp.sqrt(hilfs_Intensität_senkrecht/I0)), p0=[7,-0.37])
params_senkrecht2, cov_senkrecht2 = curve_fit(E_r_orth, Winkel[:-3], unp.nominal_values(unp.sqrt(hilfs_Intensität_senkrecht/I0)[:-3]), p0=[7,-0.37])

x_alpha = np.linspace(Winkel[0], Winkel[-1],10000)
x_alpha2 = np.linspace(Winkel[0], Winkel[-3],10000)


fig2, ax2 = plt.subplot_mosaic([["a", "a"], ["b", "b"]], sharex=True, sharey=True ,layout="constrained")

ax2["a"].errorbar(Winkel,  unp.nominal_values(unp.sqrt(hilfs_Intensität_senkrecht/I0)),xerr=0.5, yerr=unp.std_devs(unp.sqrt(hilfs_Intensität_senkrecht/I0)) , fmt="x", label="Messdaten senkrecht polarisiertes Licht")
#ax2["a"].plot(x_alpha, E_r_orth(x_alpha, 2, 0.5), label="geratener Plot")
ax2["a"].plot(x_alpha, E_r_orth(x_alpha, *params_senkrecht), "r", label ="Fit aller Werte")
ax2["a"].plot(x_alpha2, E_r_orth(x_alpha2, *params_senkrecht2), "g", label ="Fit korrigierter Werte")


params_parallel, cov_parallel = curve_fit(E_r_parallel, Winkel, unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)), p0=[7,1])
params_parallel2, cov_parallel2 = curve_fit(E_r_parallel, Winkel[:-3], unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)[:-3]), p0=[7,1])

ax2["b"].errorbar(Winkel,  unp.nominal_values(unp.sqrt(hilfs_Intensität_parallel/I0)), xerr=0.5, yerr=unp.std_devs(unp.sqrt(hilfs_Intensität_parallel/I0))  , fmt="x", label="Messdaten parallel polarisiertes Licht")
#ax2["b"].plot(x_alpha, E_r_parallel(x_alpha, 5.5,1), label="geratener Plot")
ax2["b"].plot(x_alpha, E_r_parallel(x_alpha, *params_parallel), "r", label ="Fit aller Werte")
ax2["b"].plot(x_alpha2, E_r_parallel(x_alpha2, *params_parallel2), "g", label ="Fit korrigierter Werte")

print("Das sind die shapes", x_alpha.shape, E_r_parallel(x_alpha2, *params_parallel2).shape)
ax2["a"].set(
    xlabel="Winkel/°",
    ylabel=r"$\sqrt{\frac{ I_{r_\bot} }{I_0}}$",
    )
ax2["b"].set(
    xlabel="Winkel/°",
    ylabel=r"$\sqrt{\frac{ I_{r_\parallel} }{I_0}}$",
    )
ax2["a"].legend()
ax2["b"].legend(loc="upper left")

fig2.savefig("build/I0_fit.pdf")

print("senkrecht, n, a, b:", *params_senkrecht2)
print("parallel, n, a, b:", *params_parallel2)


#Berechnung der Brechungsindizes
err_parallel  = np.sqrt(np.diag(cov_parallel2 ))
err_senkrecht = np.sqrt(np.diag(cov_senkrecht2))

n_parallel  = ufloat(params_parallel2[0] , err_parallel[0] )
n_senkrecht = ufloat(params_senkrecht2[0], err_senkrecht[0])
a_parallel  = ufloat(params_parallel2[1] , err_parallel[1] )
a_senkrecht = ufloat(params_senkrecht2[1], err_senkrecht[1])


print("Brechungsindex bei parallelem Licht:", n_parallel)
print("Brechungsindex bei senkrechtem Licht:", n_senkrecht)
print("Korrekturfaktor bei parallelem Licht:", a_parallel)
print("Korrekturfaktor bei senkrechtem Licht:", a_senkrecht)



### Brechungsindex über Gleichung 22 berechnen
n_tan = unp.tan(2*np.pi*Brewster_Winkel/360)
print("Brechungsindex über Brewsterwinkel bestimmt", n_tan)



### Abweichungen vom Literaturwert
n_lit = 3.88
delta_parallel  = 100 * abs(n_lit - n_parallel)   / n_lit
delta_senkrecht = 100 * abs(n_lit - n_senkrecht)  / n_lit
delta_tan       = 100 * abs(n_lit - n_tan)        / n_lit

print("Abweichung Brechungsindex parallel vom Literaturwert ist ", delta_parallel, "%")
print("Abweichung Brechungsindex senkrecht vom Literaturwert ist ", delta_senkrecht, "%")
print("Abweichung Brechungsindex tan vom Literaturwert ist ", delta_tan, "%")


### Berechnung des Brewsterwinkels mit Literaturwert
alpha_lit = 360 * np.arctan(n_lit) / (2*np.pi)
print("Der berechnete Literaturwert für den Brewsterwinkel beträgt ", alpha_lit, "°")

delta_alpha = 100 * abs(alpha_lit - Brewster_Winkel)/(alpha_lit)
print("Abweichung des Brewsterwinkels zum Literaturwert: ", delta_alpha, "%")

print("Intensitäten senkrecht:\n", Intensität_senkrecht)
print("Fehler senkrecht:\n", Fehler_senkrecht)
print("Intensitäten parallel: \n", Intensität_parallel)
print("Fehler parallel: \n", Fehler_parallel)