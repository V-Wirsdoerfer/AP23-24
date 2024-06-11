import matplotlib.pyplot as plt
import numpy as np

### Daten generieren

I_Anode, U = np.genfromtxt("Kennlinie2.4.txt", unpack=True)



###Testplot
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(I_Anode, U, "x", label="Testplot")
plt.show()