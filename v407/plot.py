import matplotlib.pyplot as plt
import numpy as np

#Daten generieren
Winkel, Intensität_senkrecht, drittel_senkrecht, Intensität_parallel, drittel_parallel, Fehler_senkrecht, Fehler_parallel = np.genfromtxt("Messwerte.txt", unpack=True)  #in °, A, A, A, A

#Manche Werte durch 3 teilen
drittel_senkrecht = np.asarray(drittel_senkrecht, dtype=bool)
drittel_parallel  = np.asarray(drittel_parallel,  dtype=bool)

Intensität_senkrecht[drittel_senkrecht] /= 3
Intensität_parallel[drittel_parallel]   /= 3


print("Intensität parallel:\n", np.round(Intensität_senkrecht, 7))

