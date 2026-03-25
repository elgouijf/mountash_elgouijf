import numpy as np
import matplotlib.pyplot as plt

# charger les données
data = np.loadtxt("lj.txt")

r = data[:, 0]
U = data[:, 1]

plt.figure()
plt.plot(r, U, label="Lennard-Jones")
plt.axhline(0, linestyle='--')

plt.xlabel("r")
plt.ylabel("U(r)")
plt.title("Potentiel de Lennard-Jones")
plt.legend()
plt.grid()

plt.show()