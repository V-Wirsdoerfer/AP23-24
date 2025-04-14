import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const

c = const.c
e = const.e
pi= const.pi
epsilon_0 = const.epsilon_0

#Daten generieren
x, B = np.genfromtxt("content/Magnetfeld.txt", unpack=True)




#Daten plotten
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(x,B, label="Magnetfeldstärke")
ax1.legend()
ax1.set(
    xlabel="Abstand zur Mitte in mm",
    ylabel="Magnetfeldstärke in mT"
)
fig1.savefig("Magnetfeld.pdf")


fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(x,B, label="Winkel, 1_2e18")
ax2.legend()
ax2.set(
    xlabel="Wellenlänge des Lichts in µm",
    ylabel="Winkeländerung"
)
fig2.savefig("1_2e18.pdf")


