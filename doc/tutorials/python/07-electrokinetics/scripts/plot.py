import numpy as np
import matplotlib.pyplot as plt

data_anal = np.loadtxt("eof_analytical.dat")
data_ek = np.loadtxt("eof_electrokinetics.dat")

fig1 = plt.figure()
ax = fig1.add_subplot(221)
ax.plot(data_ek[:,0], data_ek[:,1], color="r", label="electrokinetics")
ax.plot(data_anal[:,0], data_anal[:,1], color="b", label="analytical")
ax.set_xlabel("z-position")
ax.set_ylabel("density")
ax.legend(loc="best")

ax = fig1.add_subplot(222)
ax.plot(data_ek[:,0], data_ek[:,2], color="r", label="electrokinetics")
ax.plot(data_anal[:,0], data_anal[:,2], color="b", label="analytical")
ax.set_xlabel("z-position")
ax.set_ylabel("velocity")
ax.legend(loc="best")

ax = fig1.add_subplot(223)
ax.plot(data_ek[:,0], data_ek[:,3], color="r", label="electrokinetics")
ax.plot(data_anal[:,0], data_anal[:,3], color="b", label="analytical")
ax.set_xlabel("z-position")
ax.set_ylabel("shear stress xz")
ax.legend(loc="best")

plt.show()
