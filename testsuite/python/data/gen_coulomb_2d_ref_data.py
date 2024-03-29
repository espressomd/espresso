# Copyright (C) 2010-2022 The ESPResSo project
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
import scipy

box_l = 10

data = scipy.genfromtxt("coulomb_mixed_periodicity_system.data")
n = data.shape[0]

pos = data[:, 1:4]
q = data[:, 4]
forces = scipy.zeros((n, 3))
energy = 0.

q1q2 = scipy.outer(q, q)
images = 2000

for i in range(n):
    for x in range(-images, images + 1, 1):
        for y in range(-images, images + 1, 1):
            if x**2 + y**2 > images**2:
                continue
            pos_diff = pos[i] - (pos + scipy.array((x, y, 0)) * box_l)
            r = scipy.sqrt(scipy.sum(pos_diff**2, 1))
            r3 = r**3
            qq = q1q2[i, :]

            tmp = qq / r
            tmp[abs(tmp) == scipy.inf] = 0
            energy += scipy.sum(tmp)
            pref = qq / r**3
            pref = pref.reshape((n, 1))
            tmp = pos_diff * scipy.hstack((pref, pref, pref))

            forces += scipy.nan_to_num(tmp)


ids = scipy.arange(n)


forces *= -1
scipy.savetxt(
    "coulomb_mixed_periodicity_system.data",
    scipy.hstack((ids.reshape((n, 1)), pos, q.reshape((n, 1)), forces)))
