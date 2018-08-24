#ifndef OBSERVABLES_LBOBSERVABLE_HPP
#define OBSERVABLES_LBOBSERVABLE_HPP

#include <cmath>

#include "Observable.hpp"

namespace Observables {

class LBObservable : virtual public Observable {
public:
  double sampling_delta_x = 1.0, sampling_delta_y = 1.0, sampling_delta_z = 1.0;
  double sampling_offset_x = 0.0, sampling_offset_y = 0.0,
         sampling_offset_z = 0.0;
  bool allow_empty_bins = false;
  void calculate_sample_positions() {
    m_sample_positions.clear();
    if (sampling_delta_x == 0 or sampling_delta_y == 0 or sampling_delta_z == 0)
      throw std::runtime_error("Parameter delta_x/y/z must not be zero!");
    const auto n_samples_x = static_cast<size_t>(
         std::rint((box_l[0] - sampling_offset_x) / sampling_delta_x));
    const auto n_samples_y = static_cast<size_t>(
         std::rint((box_l[1] - sampling_offset_y) / sampling_delta_y));
    const auto n_samples_z = static_cast<size_t>(
         std::rint((box_l[2] - sampling_offset_z) / sampling_delta_z));
    for (size_t x = 0; x < n_samples_x; ++x) {
      for (size_t y = 0; y < n_samples_y; ++y) {
        for (size_t z = 0; z < n_samples_z; ++z) {
          m_sample_positions.push_back(sampling_offset_x + x * sampling_delta_x);
          m_sample_positions.push_back(sampling_offset_y + y * sampling_delta_y);
          m_sample_positions.push_back(sampling_offset_z + z * sampling_delta_z);
        }
      }
  }
  }
  std::vector<double> m_sample_positions;

};

}

#endif
