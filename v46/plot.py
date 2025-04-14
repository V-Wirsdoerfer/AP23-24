import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

c = const.c
e = const.e
pi= const.pi
epsilon_0 = const.epsilon_0

N_1 = 1.2e18
N_2 = 2.8e18 

#Daten generieren
x, B = np.genfromtxt("content/Magnetfeld.txt", unpack=True)
Wellenlänge, deg_0_rr, min_0_rr, deg_0_rb, min_0_rb = np.genfromtxt("content/hochrein.txt", unpack=True)    #reines GaAs ohne n-Dotierung
Wellenlänge, deg_1_rr, min_1_rr, deg_1_rb, min_1_rb = np.genfromtxt("content/1_2e18.txt", unpack=True)      #leicht n-dotiertes GaAs
Wellenlänge, deg_2_rr, min_2_rr, deg_2_rb, min_2_rb = np.genfromtxt("content/2_8e18.txt", unpack=True)      #stark n-dotiertes GaAs


#Winkeldurchschnitt bestimmen
deg_0 = 0.5*(deg_0_rr + deg_0_rb)
deg_1 = 0.5*(deg_1_rr + deg_1_rb)
deg_2 = 0.5*(deg_2_rr + deg_2_rb)


#Daten plotten
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(x,B,".", label="Magnetfeldstärke")
ax1.legend()
ax1.set(
    xlabel="Abstand zur Mitte in mm",
    ylabel="Magnetfeldstärke in mT"
)
fig1.savefig("Magnetfeld.pdf")


fig2, ax2 = plt.subplots(layout="constrained")
#ax2.plot(rad_1, label="Winkel, 1_2e18")
ax2.legend()
ax2.set(
    xlabel="Wellenlänge des Lichts in µm",
    ylabel="Winkeländerung"
)
fig2.savefig("1_2e18.pdf")


