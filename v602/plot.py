import matplotlib.pyplot as plt
import numpy as np



### Daten generieren

theta_Bragg, Rate_Bragg = np.genfromtxt("BraggBed.txt", unpack=True)
theta_Cu, Rate_Cu = np.genfromtxt("EmissionCu.txt", unpack=True)
theta_Zn, Rate_Zn = np.genfromtxt("EmissionZn.txt", unpack=True)
theta_Ga, Rate_Ga = np.genfromtxt("EmissionGa.txt", unpack=True)
theta_Br, Rate_Br = np.genfromtxt("EmissionBr.txt", unpack=True)
theta_Sr, Rate_Sr = np.genfromtxt("EmissionSr.txt", unpack=True)
theta_Zr, Rate_Zr = np.genfromtxt("EmissionZr.txt", unpack=True)


#2 theta in theta umrechnen
theta_Cu /= 2
theta_Zn /= 2
theta_Ga /= 2
theta_Br /= 2
theta_Sr /= 2
theta_Zr /= 2


print(theta_Bragg, Rate_Bragg)


### Daten plotten

# Bragbedingung
fig1, ax1 = plt.subplots(layout="constrained")
ax1.plot(theta_Bragg, Rate_Bragg, "x", label="Braggbedingung")

ax1.set(
    xlabel=r"$2 \cdot \theta/°$",
    ylabel="Zählrate",
)
ax1.legend()
fig1.savefig("build/BraggBed.pdf")

#Emissionsspektrum Cu-Röntgenröhre
fig2, ax2 = plt.subplots(layout="constrained")
ax2.plot(theta_Cu, Rate_Cu, "x", label="Emissionsspektrum Cu")

ax2.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax2.legend()
fig2.savefig("build/EmissionCu.pdf")

#Emissionsspektrum Zn Probe
fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(theta_Zn, Rate_Zn, "x", label="Emissionsspektrum Zn")

ax3.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax3.legend()
fig3.savefig("build/EmissionZn.pdf")

#Emissionsspektrum Ga Probe
fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(theta_Ga, Rate_Ga, "x", label="Emissionsspektrum Ga")

ax3.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax3.legend()
fig3.savefig("build/EmissionGa.pdf")

#Emissionsspektrum Ga Probe
fig3, ax3 = plt.subplots(layout="constrained")
ax3.plot(theta_Ga, Rate_Ga, "x", label="Emissionsspektrum Ga")

ax3.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax3.legend()
fig3.savefig("build/EmissionGa.pdf")

#Emissionsspektrum Br Probe
fig4, ax4 = plt.subplots(layout="constrained")
ax4.plot(theta_Br, Rate_Br, "x", label="Emissionsspektrum Br")

ax4.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax4.legend()
fig4.savefig("build/EmissionBr.pdf")

#Emissionsspektrum Sr Probe
fig5, ax5 = plt.subplots(layout="constrained")
ax5.plot(theta_Sr, Rate_Sr, "x", label="Emissionsspektrum Sr")

ax5.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax5.legend()
fig5.savefig("build/EmissionSr.pdf")

#Emissionsspektrum Zr Probe
fig6, ax6 = plt.subplots(layout="constrained")
ax6.plot(theta_Zr, Rate_Zr, "x", label="Emissionsspektrum Zr")

ax6.set(
    xlabel=r"$\theta/°$",
    ylabel="Zählrate",
)
ax6.legend()
fig6.savefig("build/EmissionZr.pdf")







