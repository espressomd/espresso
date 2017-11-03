#ifndef OBSERVABLES_RADIALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_RADIALFLUXDENSITYPROFILE_HPP

#include "utils.hpp"
#include "RadialProfileObservable.hpp"

namespace Observables {
class RadialFluxDensityProfile : public RadialProfileObservable {
public:
  virtual int actual_calculate(PartCfg & partCfg) override {
    double bin_volume;
    int r_bin, phi_bin, z_bin;

    for (int id : ids()) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      auto const ppos_cyl = transform_to_cylinder_coordinates(ppos, center);
      r_bin = std::floor((ppos_cyl[0] - min_r) / r_bin_size);
      phi_bin = std::floor((ppos_cyl[1] - min_phi) / phi_bin_size);
      z_bin = std::floor((ppos_cyl[2] - min_z) / z_bin_size);
      bin_volume = PI * ((min_r + (r_bin + 1) * r_bin_size) * (min_r + (r_bin + 1) * r_bin_size) - 
                         (min_r +  r_bin * r_bin_size) * (min_r + r_bin * r_bin_size)) * 
                   z_bin_size * phi_bin_size / (2 * PI);
      for (int dim = 0; dim < 3; dim++)
        last_value[3 * (r_bin * n_phi_bins * n_z_bins + phi_bin * n_z_bins + z_bin) + dim] +=
            partCfg[id].m.v[dim] / bin_volume;
    }
    return 0;
  }
  virtual int n_values() const override { return 3 * n_r_bins * n_phi_bins * n_z_bins; }
};

} // Namespace Observables

#endif
