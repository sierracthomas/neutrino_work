import matplotlib.pyplot as plt
from numpy import genfromtxt

data = genfromtxt("rp018040.tot")


plt.scatter(data.T[0], data.T[1])
plt.title("Cross section of 40Ar remaining after incident gamma at specified energies")
plt.xlabel("Energies (MeV)")
plt.ylabel("Cross section")
plt.savefig("xs_ar40.png")
plt.show()

data = genfromtxt("binary.tot")
plt.scatter(data.T[0], data.T[1])
plt.title("Binary")
plt.xlabel("Energies (MeV)")
plt.ylabel("Cross section")
plt.savefig("xs_ar40_binary.png")
plt.show()
