import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp
from scipy.stats import sem


### Konstanten und Grundgrößen ###

g = 9.81            #in m/s^2
rho_Oel = 886       #in kg/m^3
rho_L = 1.1839      #in kg/m^3 bei T=25°C
x = 0.5e-3          #in m
U = 203             #in V
d = ufloat(7.6250, 0.0051) # in mm
d *= 1e-3           #von mm in m
B = 822.599e-5      #in Pa m
p = 101325          #in Pa (Umgebungsluftdruck)
e0 = 1.602e-19      #Elementarladung in C

E=U/d


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
idx = [1, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]

auf_nom_ = []
ab_nom_ = []
auf_err_ = []
ab_err_ = []
T = []

for i in idx:
    auf, ab, temp = np.genfromtxt(f"./data/drop{i}.txt", unpack=True)
    
    auf_nom_.append(np.mean(auf))
    auf_err_.append(sem(auf))
    ab_nom_.append(np.mean(ab))
    ab_err_.append(sem(ab))

    T.append(temp[0])

auf = unp.uarray(auf_nom_, auf_err_)
ab = unp.uarray(ab_nom_, ab_err_)
T = np.asarray(T)


### in Kelvin umrechnen
T += 273.15

### 8. Tröpfchen mit Länge 1,5mm gemessen, jetzt wie anderen auf 0,5mm "normieren" ###
auf[4] /= 3
ab[4]  /= 3



### Funktionen und Formeln ###

def func_eta_Luft(T):
    return (eta_Steigung * T + eta_0)*1e-5

def Ladung(T, v_ab, v_auf):
    #if ((v_ab - v_auf)<0): 
    #    print(f"Ladung: Bei dem {i}-ten Tröpfchen ist die aufwärtsgeschwindigkeit größer als die abwärtsgeschwindigkeit.")
    q = 3 * np.pi * func_eta_Luft(T) * (v_ab + v_auf) * unp.sqrt( abs( (9 * func_eta_Luft(T) * (v_ab - v_auf) ) / (2 * g * (rho_Oel - rho_L)) ) ) / E
    return q

def Radius(T, v_ab, v_auf):
    #if ((v_ab - v_auf)<0): 
    #    print(f"Radius: Bei dem {i}-ten Tröpfchen ist die aufwärtsgeschwindigkeit größer als die abwärtsgeschwindigkeit.")
    r = 1 * unp.sqrt( abs( (9 * func_eta_Luft(T) * (v_ab - v_auf) ) / (2 * g * (rho_Oel - rho_L)) ) ) # Ganz verrückt! Das 1 * muss sein, da r sonst einen falschen arraytypen hat es wäre sonst ein Objekt
    return r

def func_q_korrigiert(q, r):
    q_korrigiert = q*(unp.pow( (1 + B/(p*r) ), (-3/2) ))
    #print("korrekturterm: ", (unp.pow( (1 + B/(p*r) ), (-3/2) )))
    return q_korrigiert

def Achsenabschnit(params, cov):
    err = np.sqrt(np.diag(cov))
    return ufloat(params[1], err[1])

### Geschwindigkeiten bestimmen
v_auf= x/auf
v_ab = x/ab



### Ladung und Radius bestimmen

q = unp.uarray(np.ones(12), np.ones(12))

for i in np.arange(12):
    q[i] = Ladung(T[i], v_ab[i], v_auf[i])
#print("q: ", q)


r = unp.uarray(np.ones(12), np.ones(12))

for i in np.arange(12):
    r[i]= Radius(T[i], v_ab[i], v_auf[i])
#print("r: ", r)


### korrigierte Ladung

q_korrigiert = unp.uarray(np.ones(12), np.ones(12))

for i in np.arange(12):
    q_korrigiert[i] = func_q_korrigiert(q[i], r[i])
#print("korrigiertes q: ", q_korrigiert)



### Versuch des plottens

fig1, ax = plt.subplots()
errx = unp.std_devs(q)
errx = unp.std_devs(q_korrigiert)
plt.errorbar(unp.nominal_values(q) + errx, np.arange(12), xerr =  errx, fmt="x" , label="Ladung unbereinigt")
plt.errorbar(unp.nominal_values(q_korrigiert) + errx, np.arange(12)+0.2, xerr =  errx, fmt="x" , label="Ladung bereinigt")

plt.axvline(x=  e0)
plt.axvline(x=2*e0)
plt.axvline(x=3*e0)
plt.axvline(x=4*e0)
ax.set(
    xlabel="Ladung in C"
)
ax.legend()
fig1.savefig("build/Ladungsauftragung.pdf")


fig2, ax = plt.subplots()
errx = unp.std_devs(q[0:4])
errx = unp.std_devs(q_korrigiert[0:4])
plt.errorbar(unp.nominal_values(q[0:4]) + errx, np.arange(4), xerr =  errx, fmt="rx" , label="Ladung unbereinigt")
plt.errorbar(unp.nominal_values(q_korrigiert[0:4]) + errx, np.arange(4)+0.2, xerr =  errx, fmt="bx" , label="Ladung bereinigt")

errx = unp.std_devs(q[5::])
errx = unp.std_devs(q_korrigiert[5::])
plt.errorbar(unp.nominal_values(q[5::]) + errx, np.arange(12)[5::], xerr =  errx, fmt="rx")
plt.errorbar(unp.nominal_values(q_korrigiert[5::]) + errx, np.arange(12)[5::]+0.2, xerr =  errx, fmt="bx")
plt.axvline(x=  e0)
plt.axvline(x=2*e0)
plt.axvline(x=3*e0)
plt.axvline(x=4*e0)
ax.set(
    xlabel="Ladung in C"
)
ax.legend()
fig2.savefig("build/Ladung_exkludiert.pdf")



### Ladung durch natürliche Zahl so lange teilen, bis

counter_korrigiert = np.zeros(12)
counter_unkorrigiert = np.zeros(12)
q_Wert_unkorrigiert = unp.uarray(np.ones(12), np.ones(12))
q_Wert_korrigiert = unp.uarray(np.ones(12), np.ones(12))

for i in np.arange(12):
    if q[i] > e0: 
        while q[i]/(counter_unkorrigiert[i]+1) > e0:
            counter_unkorrigiert[i] += 1 
            #print("counter: ", counter)
        q_Wert_unkorrigiert[i] = q[i]/counter_unkorrigiert[i]
    else: 
        print(f"Ladung des Tröpfchens Nr. {i+1} ist kleiner als e0")
        q_Wert_unkorrigiert[i] = q[i]
    
    if q_korrigiert[i] > e0: 
        while q_korrigiert[i]/(counter_korrigiert[i]+1) > e0:
            counter_korrigiert[i] += 1 
        q_Wert_korrigiert[i] = q_korrigiert[i]/counter_korrigiert[i]
            #print("counter: ", counter_korrigiert)
    else: 
        print(f"Ladung des Tröpfchens Nr. {i+1} ist kleiner als e0")
        q_Wert_korrigiert[i]= q_korrigiert[i]



### Versuch Mittelwert über lineare Regression zu berechnen ###

fig3, ax3 = plt.subplots()
plt.errorbar(np.arange(12)  , unp.nominal_values(q_Wert_unkorrigiert), yerr=unp.std_devs(q_Wert_unkorrigiert), fmt="x" , label="unkorrigierte Elementarladung")
plt.errorbar(np.arange(12)+0.2, unp.nominal_values(q_Wert_korrigiert), yerr=unp.std_devs(q_Wert_korrigiert  ), fmt="x" , label="korrigierte   Elementarladung")

err_unkorrigiert = np.asarray(unp.std_devs(q_Wert_unkorrigiert), dtype=np.float64) #polyfit kennt sonst Datentyp nicht
params_unkorrigiert, cov_unkorrigiert = np.polyfit(np.arange(12), unp.nominal_values(q_Wert_unkorrigiert), 1, w=[1/ty for ty in err_unkorrigiert], full=False, cov=True)
ax3.plot(np.arange(13), params_unkorrigiert[0]*np.ones(13)+ params_unkorrigiert[1],label="lineare Regression, unkorrigierter Werte")

err_korrigiert = np.asarray(unp.std_devs(q_Wert_korrigiert), dtype=np.float64) #polyfit kennt sonst Datentyp nicht
params_korrigiert, cov_korrigiert = np.polyfit(np.arange(12), unp.nominal_values(q_Wert_korrigiert), 1, w=[1/ty for ty in err_korrigiert], full=False, cov=True)
ax3.plot(np.arange(13), params_korrigiert[0]*np.ones(13)+ params_korrigiert[1],label="lineare Regression korrigierter Werte")

ax3.set(
    xlabel="interne Nummer des Tröpfchens",
    ylabel="Ladung / C",
)
ax3.legend()
fig3.savefig("build/e0_mean.pdf")


#rangezoomt plotten
fig4, ax4 = plt.subplots()
mask_q = [0,1,2,3,4,5,6,7,8,10]
plt.errorbar(np.arange(10)  , unp.nominal_values(q_Wert_unkorrigiert[mask_q]), yerr=unp.std_devs(q_Wert_unkorrigiert[mask_q]), fmt="x" , label="unkorrigierte Elementarladung")
plt.errorbar(np.arange(10)+0.2, unp.nominal_values(q_Wert_korrigiert[mask_q]), yerr=unp.std_devs(q_Wert_korrigiert[mask_q]  ), fmt="x" , label="korrigierte   Elementarladung")
ax4.plot(np.arange(11), params_unkorrigiert[0]*np.ones(11)+ params_unkorrigiert[1],label="lineare Regression, unkorrigierter Werte")
ax4.plot(np.arange(11), params_korrigiert[0]*np.ones(11)+ params_korrigiert[1],label="lineare Regression korrigierter Werte")
ax4.set(
    xlabel="interne Nummer des Tröpfchens",
    ylabel="Ladung / C",
)
ax4.legend()
fig4.savefig("build/e0_reduziert.pdf")



### Ladung aus params berechnen
print("Versuch mit linregess: q unkorrigiert: ", Achsenabschnit(params_unkorrigiert, cov_unkorrigiert))

print("Versuch mit linregess: q korrigiert: ", Achsenabschnit(params_korrigiert, cov_korrigiert))