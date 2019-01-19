#
# Copyright (C) 2019 The ESPResSo project
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
#

import matplotlib.pyplot as plt
plt.ion()

# get the simulated data points
sim = np.array(node_v_list)
# x-axis: in ESPResSo, LB nodes are shifted by 0.5 agrid
pos = [(i + 0.5) * agrid for i in range(len(sim))]
# analytical curve: the box width is not box_l, but box_l - 2 * wall_offset
# also, the velocity is zero beyond the walls
ana = np.array([max(0, force_density / 2. / visc
                    * ((box_l - 2 * wall_offset)**2 / 4. - (x - box_l / 2.)**2)) for x in pos])
# fit the simulated curve to the analytical curve with a least squares fit
fit = ana / sim
scaling = np.mean(fit[np.nonzero(np.invert(np.isnan(fit)))])

# plot
plt.figure(figsize=(10, 6), dpi=80)
plt.plot(pos, ana, label='analytical')
plt.plot(pos, scaling * sim, '+', label='simulated')
plt.xlabel('$x$-axis ($\AA$)', fontsize=20)
plt.ylabel('LB fluid velocity $u_y(x)$ ($\AA/s)$', fontsize=20)
plt.legend(fontsize=20)
plt.show()
