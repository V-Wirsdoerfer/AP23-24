import matplotlib.pyplot as plt
import numpy as np

x, F = np.genfromtxt("/content/Daten.txt", unpack = "True")

fig, ax = plt.subplots()

ax.plot(x,F)

plt.show()