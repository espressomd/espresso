#ifndef CORE_UNIT_TESTS_COMMON_GAUSSIAN_HPP
#define CORE_UNIT_TESTS_COMMON_GAUSSIAN_HPP

#include "Vector.hpp"

#include <boost/multi_array.hpp>

#include <cmath>

inline double gaussian(Vector3d x, Vector3d x0, double sigma) {
  return std::exp(-((x - x0).norm2() / (2. * sigma * sigma)));
}

inline Vector3d del_gaussian(Vector3d x, Vector3d x0, double sigma) {
  return -(x - x0) * gaussian(x, x0, sigma) / (sigma * sigma);
}

inline boost::multi_array<double, 3> gaussian_field(int size, Vector3d h,
                                                    Vector3d origin,
                                                    Vector3d x0, double sigma) {
  boost::multi_array<double, 3> data(Vector<3, int>{10, 10, 10});
  for (int i = 0; i < 10; i++)
    for (int j = 0; j < 10; j++)
      for (int k = 0; k < 10; k++) {
        auto const x = origin + Vector3d{i * h[0], j * h[1], k * h[2]};
        data[i][j][k] = gaussian(x, x0, sigma);
      }

  return data;
}

#endif
