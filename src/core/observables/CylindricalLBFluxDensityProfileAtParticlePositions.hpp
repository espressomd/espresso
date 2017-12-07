/*
  Copyright (C) 2016,2017 The ESPResSo project

  This file is part of ESPResSo.


  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP
#define OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP

#include "CylindricalProfileObservable.hpp"
#include "partCfg_global.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalLBFluxDensityProfileAtParticlePositions
    : public CylindricalProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override;
  virtual int n_values() const override {
    return 3 * n_r_bins * n_phi_bins * n_z_bins;
  }

private:
  virtual void do_write() override {
    // We override the implementation to actually write positions not plain
    // indices.
    static const int len_dims[4] = {n_r_bins, n_phi_bins, n_z_bins, 3};
    static const int n_dims = 4;
    static const std::array<double, 3> bin_sizes = {
        {r_bin_size(), phi_bin_size(), z_bin_size()}};
    std::array<double, 3> position;
    int index;
    int unravelled_index[4];
    std::vector<double> tmp = operator()(partCfg());
    for (auto it = tmp.begin(); it != tmp.end(); it += 3) {
      index = std::distance(tmp.begin(), it);
      ::Utils::unravel_index(len_dims, n_dims, index, unravelled_index);
      position = {
          {(static_cast<double>(unravelled_index[0]) + 0.5) * bin_sizes[0],
           (static_cast<double>(unravelled_index[1]) + 0.5) * bin_sizes[1],
           (static_cast<double>(unravelled_index[2]) + 0.5) * bin_sizes[2]}};
      m_ofile << position[0] << " " << position[1] << " " << position[2] << " "
              << *it << " " << *(it + 1) << " " << *(it + 2) << "\n";
    }
    m_ofile << std::endl;
  }
};

} // Namespace Observables

#endif
