# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import numpy as np
import matplotlib.pyplot as plt

data_anal = np.loadtxt("eof_analytical.dat")
data_ek = np.loadtxt("eof_electrokinetics.dat")

fig1 = plt.figure(figsize=(16, 4))
ax = fig1.add_subplot(131)
ax.plot(data_ek[:, 0], data_ek[:, 1], 'o',
        mfc='none', color="r", label="electrokinetics")
ax.plot(data_anal[:, 0], data_anal[:, 1], color="b", label="analytical")
ax.set_xlabel("x-position")
ax.set_ylabel("density")
ax.legend(loc="best")

ax = fig1.add_subplot(132)
ax.plot(data_ek[:, 0], data_ek[:, 2], 'o',
        mfc='none', color="r", label="electrokinetics")
ax.plot(data_anal[:, 0], data_anal[:, 2], color="b", label="analytical")
ax.set_xlabel("x-position")
ax.set_ylabel("velocity")
ax.legend(loc="best")

ax = fig1.add_subplot(133)
ax.plot(data_ek[:, 0], data_ek[:, 3], 'o',
        mfc='none', color="r", label="electrokinetics")
ax.plot(data_anal[:, 0], data_anal[:, 3], color="b", label="analytical")
ax.set_xlabel("x-position")
ax.set_ylabel("shear stress xz")
ax.legend(loc="best")

plt.show()
