import matplotlib.pyplot as plt
import numpy as np

### Literaturwerte

c_Wasser = 1497        #m/s       
c_Acryl = 2730         #m/s 
c_Aluminium = 6320     #m/s

### Daten generieren

Hoehe_Alu, Laufzeit_Alu_2MHz = np.genfromtxt("content/Aluminium2MHz.txt", unpack=True)
Hoehe_Acryl_2MHz, Laufzeit_Acryl_2MHz, U_0_Acryl_2MHz, U_Acryl_2MHz = np.genfromtxt("content/Acryl2MHz.txt", unpack=True)
Hoehe_Acryl_1MHz, U_0_Acryl_1MHz, U_Acryl_1MHz = np.genfromtxt("content/Acryl1MHz.txt", unpack=True)

### Schallgeschwindigkeit bestimmen 

### Aluminium 

fig, ax1 = plt.subplots(1, 1, layout="constrained")
ax1.plot(Laufzeit_Alu_2MHz, 2 * Hoehe_Alu, "rx", label="c_Alu")
params_WegZeitAlu, cov_WegZeitAlu = np.polyfit(Laufzeit_Alu_2MHz, 2 * Hoehe_Alu, deg=1, cov=True)
print(params_WegZeitAlu[0] * 100) 