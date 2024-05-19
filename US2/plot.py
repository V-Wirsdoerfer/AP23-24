import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 

###Grundgrößen und Konstanten
c_Acryl = 2730      #in m/s
c_Wasser = 1497     #in m/s
c_Linse = 2500      #in m/s
c_Glaskörper = 1410 #in m/s
H_Acrylblock = ufloat(0.079, 0.000025)  #in m
dm_Kegel = ufloat(45.5e-3, 0.000025)      #in m


### Daten generieren
H_top, H_bot = np.genfromtxt("./content/acrylblock.txt", unpack=True)         #in mm und mm
t_top, t_bot = np.genfromtxt("./content/Laufzeit_Acryl.txt", unpack=True)     #in µs und µs
Herz_T, Herz_Amp = np.genfromtxt("./content/Herz.txt", unpack=True)           #in  s und µs


### Daten in SI umrechnen
H_top *= 1e-3       #von mm in m
H_bot *= 1e-3       #von mm in m

t_top *= 1e-6       #von µs in s
t_bot *= 1e-6       #von µs in s

Herz_Amp *= 1e-6    #von µs in s

### Fehler einbeziehen 

H_top = unp.uarray(H_top, 0.025e-3)
H_bot = unp.uarray(H_bot, 0.025e-3)

### Anpassungsschicht berechnen

### Laufzeitmessung top

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(t_top, 2 * unp.nominal_values(H_top), "rx", label="Messdaten")
params_top, cov_top = np.polyfit(t_top, 2 * unp.nominal_values(H_top), deg=1, cov=True)
ax.plot(t_top, params_top[0] * t_top + params_top[1], label="Ausgleichsgerade")
ax.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $x$ / m",
    xlim = (0, 4.5e-5),
    ylim = (0, 0.13)
)
print("Schallgeschwindigkeit der top-Messung: ", params_top[0])
ax.legend()
fig.savefig("build/Schallgeschwindigkeit_top.pdf")

### Laufzeitmessung bottom

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(t_bot, 2 * unp.nominal_values(H_bot), "rx", label="Messdaten")
params_bot, cov_bot = np.polyfit(t_bot, 2 * unp.nominal_values(H_bot), deg=1, cov=True)
ax.plot(t_bot, params_bot[0] * t_bot + params_bot[1], label="Ausgleichsgerade")
ax.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $x$ / m",
#    xlim = (0, 4.5e-5),
#    ylim = (0, 0.13)
)
print("Schallgeschwindigkeit der bottom-Messung: ", params_bot[0])
ax.legend()
fig.savefig("build/Schallgeschwindigkeit_bottom.pdf")

Mittelwert_Schall = 0.5 * (params_top[0] + params_bot[0])
print("Mittelwert Schallgeschwindigkeit: ", Mittelwert_Schall)

### Laufzeitkorrektur

H_top_korr = 0.5 * t_top * c_Acryl 
d_Anp_top = H_top_korr - H_top
#print("Breite der Anpassungsschicht nach top: ", np.mean(d_Anp_top))

H_bot_korr = 0.5 * t_bot * c_Acryl 
d_Anp_bot = H_bot_korr - H_bot
#print("Breite der Anpassungsschicht nach bot: ", np.mean(d_Anp_bot))

d_Anp = 0.5 * (np.mean(d_Anp_bot) + np.mean(d_Anp_top))
print("Mittelwert der Anpassungsschicht: ", d_Anp)

t_Anp = d_Anp / c_Acryl
t_eff_top = t_top - 2 * t_Anp
t_eff_bot = t_bot - 2 * t_Anp
#print("Effektive Zeit im Acrylblock top: ", t_eff_top)
#print("Effektive Zeit im Acrylblock bot: ", t_eff_bot)


### Berechnung der Lochdurchmesser

d_top = c_Acryl * 0.5 * t_eff_top
d_bot = c_Acryl * 0.5 * t_eff_bot
dm_Loecher = H_Acrylblock - (d_top + d_bot)
#print("Das ist die Lage der Löcher aus der top-Perspektive: ", d_top)
#print("Das ist die Lage der Löcher aus der bottom-Perspektive: ", d_bot)
print("\nDas sind die Durchmesser der Löcher mit Schall in mm: \n", dm_Loecher*1e3)

d_Messschieber = H_Acrylblock - H_top - H_bot
print("\nDas sind die Durchmesser der Löcher mit Messschieber in mm: \n", d_Messschieber*1e3)

Abweichung = abs(dm_Loecher - d_Messschieber)/ d_Messschieber
print("\nAbweichung der Messmethoden relativ %: \n", Abweichung*100)
print("\nAbweichung der Messmethoden absolut in mm: \n" , abs(dm_Loecher - d_Messschieber)*1e3 )

### Hervolumen bestimmen

T_Herz = np.mean(Herz_T)
f_Herz = 1 / T_Herz
h_Herz = np.mean(Herz_Amp)
EDS = (1 / 3) * np.pi * (dm_Kegel / 2)**2 * h_Herz
HZV = EDS * f_Herz 
print("Wert für das Herzvolumen: ", HZV)
### Berechnung Schallgeschwindigkeit korrigiert

### Laufzeitmessung top_eff

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(unp.nominal_values(t_eff_top), 2 * unp.nominal_values(H_top), "rx", label="Messdaten")
params_eff_top, cov_eff_top = np.polyfit(unp.nominal_values(t_eff_top), 2 * unp.nominal_values(H_top), deg=1, cov=True)
ax.plot(unp.nominal_values(t_eff_top), params_eff_top[0] * unp.nominal_values(t_eff_top) + params_eff_top[1], label="Ausgleichsgerade")
ax.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $x$ / m",
    xlim = (0, 4.5e-5),
    ylim = (0, 0.13)
)
print("Schallgeschwindigkeit der top-Messung: ", params_eff_top[0])
ax.legend()
fig.savefig("build/Schallgeschwindigkeit_eff_top.pdf")

### Laufzeitmessung eff_bottom

fig, ax = plt.subplots(1, 1, layout="constrained")
ax.plot(unp.nominal_values(t_eff_bot), 2 * unp.nominal_values(H_bot), "rx", label="Messdaten")
params_eff_bot, cov_eff_bot = np.polyfit(unp.nominal_values(t_eff_bot), 2 * unp.nominal_values(H_bot), deg=1, cov=True)
ax.plot(unp.nominal_values(t_eff_bot), params_eff_bot[0] * unp.nominal_values(t_eff_bot) + params_eff_bot[1], label="Ausgleichsgerade")
ax.set(
    xlabel = r"Laufzeit $t$ / s",
    ylabel = r"Strecke $x$ / m",
#    xlim = (0, 4.5e-5),
#    ylim = (0, 0.13)
)
print("Schallgeschwindigkeit der eff_bottom-Messung: ", params_eff_bot[0])
ax.legend()
fig.savefig("build/Schallgeschwindigkeit_eff_bottom.pdf")
Schall_Mittel_korr = 0.5 * (params_eff_top[0] + params_eff_bot[0])
print(Schall_Mittel_korr)
#print(t_top - t_eff_top)