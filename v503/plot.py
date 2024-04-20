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


#auf1, ab1, T1 = np.genfromtxt("./data/drop1.txt", unpack=True)      #in s, s, °C
#auf4, ab4, T4 = np.genfromtxt("./data/drop4.txt", unpack=True)      #in s, s, °C
#auf6, ab6, T6 = np.genfromtxt("./data/drop6.txt", unpack=True)      #in s, s, °C
#auf7, ab7, T7 = np.genfromtxt("./data/drop7.txt", unpack=True)      #in s, s, °C
#auf8, ab8, T8 = np.genfromtxt("./data/drop8.txt", unpack=True)      #in s, s, °C
#auf9, ab9, T9 = np.genfromtxt("./data/drop9.txt", unpack=True)      #in s, s, °C
#auf10, ab10, T10 = np.genfromtxt("./data/drop10.txt", unpack=True)  #in s, s, °C    
#auf11, ab11, T11 = np.genfromtxt("./data/drop11.txt", unpack=True)  #in s, s, °C    
#auf12, ab12, T12 = np.genfromtxt("./data/drop12.txt", unpack=True)  #in s, s, °C    
#auf13, ab13, T13 = np.genfromtxt("./data/drop13.txt", unpack=True)  #in s, s, °C    
#auf14, ab14, T14 = np.genfromtxt("./data/drop14.txt", unpack=True)  #in s, s, °C    
#auf15, ab15, T15 = np.genfromtxt("./data/drop15.txt", unpack=True)  #in s, s, °C    
#
#### Mittelwerte berechnen von den Zeiten ###
#
#auf_nom = np.ones(12)
#auf_nom[0] = np.mean(auf1)
#auf_nom[1] = np.mean(auf4)
#auf_nom[2] = np.mean(auf6)
#auf_nom[3] = np.mean(auf7)
#auf_nom[4] = np.mean(auf8)
#auf_nom[5] = np.mean(auf9)
#auf_nom[6] = np.mean(auf10)
#auf_nom[7] = np.mean(auf11)
#auf_nom[8] = np.mean(auf12)
#auf_nom[9] = np.mean(auf13)
#auf_nom[10] = np.mean(auf14)
#auf_nom[11] = np.mean(auf15)
#
#auf_err = np.ones(12)
#auf_err[0] = sem(auf1)
#auf_err[1] = sem(auf4)
#auf_err[2] = sem(auf6)
#auf_err[3] = sem(auf7)
#auf_err[4] = sem(auf8)
#auf_err[5] = sem(auf9)
#auf_err[6] = sem(auf10)
#auf_err[7] = sem(auf11)
#auf_err[8] = sem(auf12)
#auf_err[9] = sem(auf13)
#auf_err[10] = sem(auf14)
#auf_err[11] = sem(auf15)
#
#auf = unp.uarray(auf_nom, auf_err)
#
#
#ab_nom = np.ones(12)
#ab_nom[0] = np.mean(ab1)
#ab_nom[1] = np.mean(ab4)
#ab_nom[2] = np.mean(ab6)
#ab_nom[3] = np.mean(ab7)
#ab_nom[4] = np.mean(ab8)
#ab_nom[5] = np.mean(ab9)
#ab_nom[6] = np.mean(ab10)
#ab_nom[7] = np.mean(ab11)
#ab_nom[8] = np.mean(ab12)
#ab_nom[9] = np.mean(ab13)
#ab_nom[10] = np.mean(ab14)
#ab_nom[11] = np.mean(ab15)
#
#ab_err = np.ones(12)
#ab_err[0] = sem(ab1)
#ab_err[1] = sem(ab4)
#ab_err[2] = sem(ab6)
#ab_err[3] = sem(ab7)
#ab_err[4] = sem(ab8)
#ab_err[5] = sem(ab9)
#ab_err[6] = sem(ab10)
#ab_err[7] = sem(ab11)
#ab_err[8] = sem(ab12)
#ab_err[9] = sem(ab13)
#ab_err[10] = sem(ab14)
#ab_err[11] = sem(ab15)
#
#ab = unp.uarray(ab_nom, ab_err)
#
#T = [T1[0], T4[0], T6[0], T7[0], T8[0], T9[0], T10[0], T11[0], T12[0], T13[0], T14[0], T15[0]]



### in Kelvin umrechnen
# for i in np.arange(12):
#     T[i] += 273.15
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

def my_mean(q):
    sum1 = 0
    sum2 = 0
    err = 0

    for i in np.arange(12):
        sum1 += q[i].n/(q[i].s**2)
        sum2 += 1     /(q[i].s**2)

    for i in np.arange(12):
        err += q[i].s**2

    value = ufloat(sum1/sum2, np.sqrt(err))
    return value


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


#print("richtige Werte unkorrigiert: ",q_Wert_unkorrigiert)
e0_unkorrigiert = ufloat(np.mean(unp.nominal_values(q_Wert_unkorrigiert)), sem(unp.nominal_values(q_Wert_unkorrigiert)))
print("richtige Werte unkorrigiert mean: ", e0_unkorrigiert)
print("Versuch richtigen Fehler zu bestimmen unkorrigiert: ", my_mean(q_Wert_unkorrigiert))


mask_korrigiert = np.array(counter_korrigiert, dtype="bool")
#print("richtige Werte korrigiert: ",q_Wert_korrigiert)
e0_korrigiert = ufloat(np.mean(unp.nominal_values(q_Wert_korrigiert)), sem(unp.nominal_values(q_Wert_korrigiert)))
print("richtige Werte korrigiert mean: ", e0_korrigiert)
print("Versuch richtigen Fehler zu bestimmen korrigiert: ", my_mean(q_Wert_korrigiert))
