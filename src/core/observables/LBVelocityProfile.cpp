#include "lb.hpp"
#include "utils/Histogram.hpp"
#include "LBVelocityProfile.hpp"

namespace Observables {

std::vector<double> LBVelocityProfile::
operator()(PartCfg &partCfg) const {
  std::array<size_t, 3> n_bins{{static_cast<size_t>(n_x_bins),
                                static_cast<size_t>(n_y_bins),
                                static_cast<size_t>(n_z_bins)}};
  std::array<std::pair<double, double>, 3> limits{
      {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
       std::make_pair(min_z, max_z)}};
  Utils::Histogram<double, 3> histogram(n_bins, 3, limits);
  // First collect all positions (since we want to call the LB function to
  // get the fluid velocities only once).
  std::vector<double> velocities(m_sample_positions.size());
  if (lattice_switch & LATTICE_LB_GPU) {
#if defined(LB_GPU)
    lb_lbfluid_get_interpolated_velocity_at_positions(
        m_sample_positions.data(), velocities.data(),
        m_sample_positions.size() / 3);
#endif
  } else if (lattice_switch & LATTICE_LB) {
#if defined(LB)
    for (size_t ind=0; ind < m_sample_positions.size(); ind +=3) {
      Vector3d pos_tmp = {m_sample_positions[ind + 0],
                           m_sample_positions[ind + 1],
                           m_sample_positions[ind + 2]};
      lb_lbfluid_get_interpolated_velocity(pos_tmp, &(velocities[ind + 0]));
    }
#endif
  } else {
    throw std::runtime_error("Either CPU LB or GPU LB has to be active for this observable to work.");
  }
  for (size_t ind = 0; ind < m_sample_positions.size(); ind += 3) {
    const Vector3d position = {
        {m_sample_positions[ind + 0], m_sample_positions[ind + 1], m_sample_positions[ind + 2]}};
    const Vector3d velocity = {
        {velocities[ind + 0], velocities[ind + 1], velocities[ind + 2]}};
    histogram.update(position, velocity);
  }
  auto hist_tmp = histogram.get_histogram();
  auto const tot_count = histogram.get_tot_count();
  for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
    if (tot_count[ind] == 0 and not allow_empty_bins) {
      auto const error = "Decrease sampling delta(s), bin " + std::to_string(ind) + " has no hit";
      throw std::runtime_error(error);
    }
    if (tot_count[ind] > 0) {
      hist_tmp[ind] /= tot_count[ind];
    }
  }
  return hist_tmp;
}

}
