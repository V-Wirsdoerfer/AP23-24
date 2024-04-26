from cProfile import label
from enum import auto
import matplotlib.pyplot as plt
import numpy as np


counts_hist = np.genfromtxt("content/Verteilung.txt", unpack = True)

### Zufallswerte erstellen
rng = np.random.default_rng(1)
mean = np.mean(counts_hist) #von wiki geklaut die Formel
poisson_Verteilung = rng.poisson(mean, size=100)
gauss_Verteilung = rng.normal(mean, 100, size=100)



#counts mit Histogramm plotten
fig2, axs = plt.subplots(2, 2, layout="constrained")
ax = np.ravel(axs)


#Histogramme in schleife erstellen
for i in np.arange(4):
    ax[i].hist(counts_hist, label="Messung", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].hist(poisson_Verteilung, label="Poisson", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].hist(gauss_Verteilung, label="Gauss", bins=(10 + 2*i), histtype="step", range=[1700,2050])
    ax[i].legend()
    ax[i].set_title(f"{10+2*i} - Bins")

fig2.savefig("build/Verteilung.pdf")

