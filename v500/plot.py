from turtle import color
import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp

### Daten generieren 
U_Blau_voll, I_Blau_voll = np.genfromtxt("blau_voll.txt", unpack=True)
U_Blau_halb, I_Blau_halb = np.genfromtxt("blau_halb.txt", unpack=True)
U_grün, I_grün, U_gelb, I_gelb, U_violett, I_violett = np.genfromtxt("grün_gelb_violett.txt", unpack=True)



# Daten in SI umrechnen
for I in (I_Blau_voll, I_Blau_halb, I_grün, I_gelb, I_violett):
    I *= 1e-9           #von nA nach A


#Fehler definieren
I_err_halb = 0.05e-9            #in A
U_err_halb = 0.01               #in V

I_err_voll_bunt = 0.0005*1e-9     # in A
U_err_voll_bunt = 0.005         # in V 



### Plot volle intensität.
fig1, ax1 = plt.subplots(layout="constrained")

ax1.errorbar(U_Blau_voll, I_Blau_voll, xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="bx", label="volle Intensität Messdaten")
subax1 = ax1.inset_axes([0.385, 0.11, 0.6, 0.4])
subax1.errorbar(U_Blau_voll[0:20], I_Blau_voll[0:20], xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="bx")

ax1.set(
    xlabel=r"$U$/V",
    ylabel=r"I/A",
#    xlim=[-1.3,2.2],
 #   ylim=[-0.2e-9,12e-9],
)
subax1.set(
    xlabel=r"$U$/V",
    ylabel=r"I/A",
)
ax1.legend()

fig1.savefig("build/blau_voll.pdf")


### Plot halbe Intensität
fig2, ax2 = plt.subplots(layout = "constrained")

ax2.errorbar(U_Blau_halb, I_Blau_halb, xerr=U_err_halb, yerr=I_err_halb, fmt="x", color="cyan", label="halbe Intensität Messdaten")
ax2.set(
    xlabel=r"$U$/V",
    ylabel=r"I/A",
#    xlim=[-1.3,2.2],
 #   ylim=[-0.2e-9,12e-9],
)
ax2.legend()

fig2.savefig("build/blau_halb.pdf")


### Beide Blaue Intensitäten auf einmal in Graphik von Blau halb hinein
ax2.errorbar(U_Blau_voll, I_Blau_voll, xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="x", color="blue", label="volle Intensität Messdaten")
ax2.legend()
fig2.savefig("build/blau_voll_halb.pdf")


### Drei anderen Farben
fig3, ax3 = plt.subplots(layout="constrained")

ax3.errorbar(U_grün, I_grün, xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="x", color="chartreuse", label="grüne Intensität Messdaten")
ax3.errorbar(U_gelb, I_gelb, xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="x", color="gold", label="gelbe Intensität Messdaten")
ax3.errorbar(U_violett, I_violett, xerr=U_err_voll_bunt, yerr=I_err_voll_bunt, fmt="x", color="purple", label="violette Intensität Messdaten")

ax3.legend()
ax3.set(
    xlabel=r"$U$/V",
    ylabel=r"I/A",
)

fig3.savefig("build/bunt.pdf")


### Ausgleichsgerade plotten

def Ausgleichsgerade(x, y, err):
    params, cov = np.polyfit(x,y, deg=1, cov=True, w=[1/err for i in range(len(x))])
    return params, cov

def Steigung(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[0], err[0])

def Achsenabschnitt(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[1], err[1])


params_blau, cov_blau = Ausgleichsgerade(U_Blau_voll[0:20], I_Blau_voll[0:20], I_err_voll_bunt)
params_grün, cov_grün = Ausgleichsgerade(U_grün, I_grün, I_err_voll_bunt)
params_gelb, cov_gelb = Ausgleichsgerade(U_gelb[0:-3], I_gelb[0:-3], I_err_voll_bunt)
params_violett, cov_violett = Ausgleichsgerade(U_violett, I_violett, I_err_voll_bunt)

Steigung_blau = Steigung(params_blau, cov_blau)
Steigung_grün = Steigung(params_grün, cov_grün)
Steigung_gelb = Steigung(params_gelb, cov_gelb)
Steigung_violett = Steigung(params_violett, cov_violett)

x_blau = np.linspace(U_Blau_voll[0], U_Blau_voll[20])
x_grün = np.linspace(U_grün[0], U_grün[-1])
x_gelb = np.linspace(U_gelb[0], U_gelb[-3])
x_violett = np.linspace(U_violett[0], U_violett[-1])

#blau
ax1.plot(x_blau, x_blau*Steigung_blau.n + Achsenabschnitt(params_blau, cov_blau).n, color="xkcd:robin's egg blue", label = "Ausgleichsgerade im linearen Bereich" )
subax1.plot(x_blau, x_blau*Steigung_blau.n + Achsenabschnitt(params_blau, cov_blau).n, color="xkcd:robin's egg blue", label = "Ausgleichsgerade im linearen Bereich" )
ax1.legend()
fig1.savefig("build/blau_voll_fit.pdf")

#grün
ax3.plot(x_grün, x_grün*Steigung_grün.n + Achsenabschnitt(params_grün, cov_grün).n, color="xkcd:seafoam green", label = "Ausgleichsgerade im linearen Bereich" )
#gelb
ax3.plot(x_gelb, x_gelb*Steigung_gelb.n + Achsenabschnitt(params_gelb, cov_gelb).n, color="xkcd:yellow", label = "Ausgleichsgerade im linearen Bereich" )
#violett
ax3.plot(x_violett, x_violett*Steigung_violett.n + Achsenabschnitt(params_violett, cov_violett).n, color="fuchsia", label = "Ausgleichsgerade im linearen Bereich" )

ax3.legend()
fig3.savefig("build/bunt_fit.pdf")


#Grenzspannungen berechnen
def Grenzspannung(params, cov):    #Nullstelle der Geradengleichung berechnen
    '''
    Stellt die Geradengleichung nach der Grenzspannung um
    '''
    m = Steigung(params, cov)
    b = Achsenabschnitt(params, cov)
    return (-b/m)

UG_blau = Grenzspannung(params_blau, cov_blau)
UG_grün = Grenzspannung(params_grün, cov_grün)
UG_gelb = Grenzspannung(params_gelb, cov_gelb)
UG_violett = Grenzspannung(params_violett, cov_violett)

name=["UG_blau", "UG_grün", "UG_gelb", "UG_violett"]
UG  =[ UG_blau ,  UG_grün ,  UG_gelb ,  UG_violett]

for U, n in zip(UG, name):
    print(n, "=", U, "V")



# Grenzspannungen plotten
fig4, ax4 = plt.subplots(layout="constrained")
ax4.errorbar(UG_violett.n, 404.7e-9, xerr=UG_violett.s, fmt="x", color="purple", label=r"$U_G$ violett")
ax4.errorbar(UG_blau.n, 435.8e-9, xerr=UG_blau.s, fmt="x", color="blue", label=r"$U_G$ blau")
ax4.errorbar(UG_grün.n, 546e-9, xerr=UG_grün.s, fmt="x", color="chartreuse", label=r"$U_G$ grün")
ax4.errorbar(UG_gelb.n, 577e-9, xerr=UG_gelb.s, fmt="x", color="gold", label=r"$U_G$ gelb")
#fit
UG =[UG_violett, UG_blau, UG_grün, UG_gelb]
UG_n =[UG_violett.n, UG_blau.n, UG_grün.n, UG_gelb.n] # muss extra definiert werden, sonst probleme bei polyfit
Wellenlaenge = [404.7e-9, 435.8e-9, 546e-9, 577e-9]
params_UG, cov_UG = np.polyfit(unp.nominal_values(UG), Wellenlaenge, deg=1, cov=True, w=[1/ty for ty in UG_n] )
x_UG = np.linspace(UG_n[0], UG_n[-1])
ax4.plot(x_UG, x_UG*Steigung(params_UG, cov_UG).n + Achsenabschnitt(params_UG, cov_UG).n, label="Ausgleichsgerade", color="orange")

ax4.set(
    ylabel=r"$\lambda$/m",
    xlabel=r"$U_G$/V",
)
ax4.legend()

fig4.savefig("build/Grenzspannung.pdf")

print("Das plancksche Wirkungsquantum:\nh =", Steigung(params_UG, cov_UG))



