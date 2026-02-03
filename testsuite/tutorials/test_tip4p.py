#
# Copyright (C) 2026 The ESPResSo project
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

import unittest as ut
import importlib_wrapper
import numpy as np
import scipy.signal
import espressomd


mesh_size = 60 if espressomd.gpu_available() else 48
tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/mlip-water/01_TIP4P_water.py",
    rdf_samples=80, p3m_params={"cao": 7, "mesh": 3 * [mesh_size]})


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test(self):
        def get_peaks_at(xdata, ydata, targets, comparator):
            """
            Determine all local extrema, but only return the subset of extrema
            that are near the expected values (targets). This makes extrema
            detection more robust against noise in the data.
            """
            idx = scipy.signal.argrelextrema(ydata, comparator)
            positions = np.array(len(xdata) * [np.inf])
            positions[idx] = xdata[idx]
            assert len(positions) >= len(targets)
            best_matches = []
            for target in targets:
                dist = np.abs(positions - target)
                index = np.argmin(dist)
                best_matches.append(index)
                positions[index] = np.inf
            return (np.array(best_matches),)

        targets_maxima = [2.8, 4.5, 6.75]
        targets_minima = [3.4, 5.5]
        ref_rs = tutorial.df["rs"].to_numpy()
        ref_rdf = tutorial.df["rdf"].to_numpy()
        sim_rs = 10. * tutorial.bin_centers
        sim_rdf = tutorial.rdf
        sim_rdf_smooth = scipy.signal.savgol_filter(sim_rdf, 6, 2)
        idx = get_peaks_at(ref_rs, ref_rdf, targets_maxima, np.greater)
        ref_maxima_x = ref_rs[idx]
        ref_maxima_y = ref_rdf[idx]
        idx = get_peaks_at(sim_rs, sim_rdf_smooth, targets_maxima, np.greater)
        sim_maxima_x = sim_rs[idx]
        sim_maxima_y = sim_rdf[idx]
        idx = get_peaks_at(ref_rs, ref_rdf, targets_minima, np.less)
        ref_minima_x = ref_rs[idx]
        ref_minima_y = ref_rdf[idx]
        idx = get_peaks_at(sim_rs, sim_rdf_smooth, targets_minima, np.less)
        sim_minima_x = sim_rs[idx]
        sim_minima_y = sim_rdf[idx]
        tol = {"rtol": 0.1}
        # compare peaks position and magnitude in the structural region
        np.testing.assert_allclose(sim_maxima_x[1:], ref_maxima_x[1:], **tol)
        np.testing.assert_allclose(sim_maxima_y[1:], ref_maxima_y[1:], **tol)
        np.testing.assert_allclose(sim_minima_x, ref_minima_x, **tol)
        np.testing.assert_allclose(sim_minima_y, ref_minima_y, **tol)
        # first peak should be 15% larger in magnitude with TIP4P/2005
        np.testing.assert_allclose(sim_maxima_x[:1], ref_maxima_x[:1], **tol)
        np.testing.assert_allclose(
            sim_maxima_y[:1], 1.15 * ref_maxima_y[:1], **tol)


if __name__ == "__main__":
    ut.main()
