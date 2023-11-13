import matplotlib.pyplot as plt
import numpy as np

x, F = np.genfromtxt("content/Daten.txt", unpack = True)
x_dots = np.linspace(0, 0.6, 1000)


#lineare regression
params, covariance_matrix = np.polyfit(x, F, deg=1, cov=True)

errors =np.sqrt(np.diag(covariance_matrix))

fig, ax = plt.subplots()

#Messwerte plotten
ax.plot(x,F, "x", label="Messdaten")

#Regressionsgerade plotten
ax.plot(x_dots, params[0]*x_dots + params[1], label = "lineare Regression")

#Plot schön machen
ax.set(
    xlabel = r"$\text{Auslenkung} \; \increment x$",
    ylabel = r"$\text{Kraft} \; F$",
)
ax.legend()


fig.savefig("build/plot.pdf")