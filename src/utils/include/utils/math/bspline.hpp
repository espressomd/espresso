/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef UTILS_MATH_BSPLINE_HPP
#define UTILS_MATH_BSPLINE_HPP

#include "sqr.hpp"
#include <utils/device_qualifier.hpp>

#include <stdexcept>
#include <type_traits>

namespace Utils {
template <int order, typename T>
DEVICE_QUALIFIER auto bspline(int i, T x)
    -> std::enable_if_t<(order > 0) && (order <= 7), T> {
  DEVICE_ASSERT(i < order);
  DEVICE_ASSERT(x >= -0.5);
  DEVICE_ASSERT(x <= 0.5);

  switch (order) {
  case 1:
    return 1.0;
  case 2: {
    switch (i) {
    case 0:
      return 0.5 - x;
    case 1:
      return 0.5 + x;
    }
  }
  case 3: {
    switch (i) {
    case 0:
      return 0.5 * sqr(0.5 - x);
    case 1:
      return 0.75 - sqr(x);
    case 2:
      return 0.5 * sqr(0.5 + x);
    }

  case 4: {
    switch (i) {
    case 0:
      return (1.0 + x * (-6.0 + x * (12.0 - x * 8.0))) / 48.0;
    case 1:
      return (23.0 + x * (-30.0 + x * (-12.0 + x * 24.0))) / 48.0;
    case 2:
      return (23.0 + x * (30.0 + x * (-12.0 - x * 24.0))) / 48.0;
    case 3:
      return (1.0 + x * (6.0 + x * (12.0 + x * 8.0))) / 48.0;
    }
  }
  case 5: {
    switch (i) {
    case 0:
      return (1.0 + x * (-8.0 + x * (24.0 + x * (-32.0 + x * 16.0)))) / 384.0;
    case 1:
      return (19.0 + x * (-44.0 + x * (24.0 + x * (16.0 - x * 16.0)))) / 96.0;
    case 2:
      return (115.0 + x * x * (-120.0 + x * x * 48.0)) / 192.0;
    case 3:
      return (19.0 + x * (44.0 + x * (24.0 + x * (-16.0 - x * 16.0)))) / 96.0;
    case 4:
      return (1.0 + x * (8.0 + x * (24.0 + x * (32.0 + x * 16.0)))) / 384.0;
    }
  }
  case 6: {
    switch (i) {
    case 0:
      return (1.0 +
              x * (-10.0 + x * (40.0 + x * (-80.0 + x * (80.0 - x * 32.0))))) /
             3840.0;
    case 1:
      return (237.0 +
              x * (-750.0 +
                   x * (840.0 + x * (-240.0 + x * (-240.0 + x * 160.0))))) /
             3840.0;
    case 2:
      return (841.0 +
              x * (-770.0 +
                   x * (-440.0 + x * (560.0 + x * (80.0 - x * 160.0))))) /
             1920.0;
    case 3:
      return (841.0 +
              x * (+770.0 +
                   x * (-440.0 + x * (-560.0 + x * (80.0 + x * 160.0))))) /
             1920.0;
    case 4:
      return (237.0 +
              x * (750.0 +
                   x * (840.0 + x * (240.0 + x * (-240.0 - x * 160.0))))) /
             3840.0;
    case 5:
      return (1.0 +
              x * (10.0 + x * (40.0 + x * (80.0 + x * (80.0 + x * 32.0))))) /
             3840.0;
    }
  }
  case 7: {
    switch (i) {
    case 0:
      return (1.0 +
              x * (-12.0 +
                   x * (60.0 + x * (-160.0 +
                                    x * (240.0 + x * (-192.0 + x * 64.0)))))) /
             46080.0;
    case 1:
      return (361.0 + x * (-1416.0 +
                           x * (2220.0 +
                                x * (-1600.0 +
                                     x * (240.0 + x * (384.0 - x * 192.0)))))) /
             23040.0;
    case 2:
      return (10543.0 +
              x * (-17340.0 +
                   x * (4740.0 +
                        x * (6880.0 +
                             x * (-4080.0 + x * (-960.0 + x * 960.0)))))) /
             46080.0;
    case 3:
      return (5887.0 + x * x * (-4620.0 + x * x * (1680.0 - x * x * 320.0))) /
             11520.0;
    case 4:
      return (10543.0 +
              x * (17340.0 +
                   x * (4740.0 +
                        x * (-6880.0 +
                             x * (-4080.0 + x * (960.0 + x * 960.0)))))) /
             46080.0;
    case 5:
      return (361.0 +
              x * (1416.0 +
                   x * (2220.0 +
                        x * (1600.0 +
                             x * (240.0 + x * (-384.0 - x * 192.0)))))) /
             23040.0;
    case 6:
      return (1.0 +
              x * (12.0 +
                   x * (60.0 +
                        x * (160.0 + x * (240.0 + x * (192.0 + x * 64.0)))))) /
             46080.0;
    }
  }
  }
  }

  DEVICE_THROW(std::runtime_error("Internal interpolation error."));
  return T{};
}

template <class T> auto bspline(int i, T x, int k) {
  switch (k) {
  case 1:
    return bspline<1>(i, x);
  case 2:
    return bspline<2>(i, x);
  case 3:
    return bspline<3>(i, x);
  case 4:
    return bspline<4>(i, x);
  case 5:
    return bspline<5>(i, x);
  case 6:
    return bspline<6>(i, x);
  case 7:
    return bspline<7>(i, x);
  }

  return 0.0;
}

template <int order, typename T = double> inline T bspline_d(int i, T x) {
  static_assert(order <= 7, "");
  DEVICE_ASSERT(i < order);
  DEVICE_ASSERT(x >= -0.5);
  DEVICE_ASSERT(x <= 0.5);

  switch (order - 1) {
  case 0:
    return 0.;
  case 1:
    switch (i) {
    case 0:
      return -1.0;
    case 1:
      return 1.0;
    }
  case 2:
    switch (i) {
    case 0:
      return x - 0.5;
    case 1:
      return -2.0 * x;
    case 2:
      return x + 0.5;
    }
  case 3:
    switch (i) {
    case 0:
      return (-1.0 + x * (4.0 + x * (-4.0))) / 8.0;
    case 1:
      return (-5.0 + x * (-4.0 + x * (12.0))) / 8.0;
    case 2:
      return (5.0 + x * (-4.0 + x * (-12.0))) / 8.0;
    case 3:
      return (1.0 + x * (4.0 + x * (4.0))) / 8.0;
    }
  case 4:
    switch (i) {
    case 0:
      return (-1.0 + x * (6.0 + x * (-12.0 + x * (8.0)))) / 48.0;
    case 1:
      return (-11.0 + x * (12.0 + x * (12.0 + x * (-16.0)))) / 24.0;
    case 2:
      return (x * (-5.0 + x * x * 4.0)) / 4.0;
    case 3:
      return (11.0 + x * (12.0 + x * (-12.0 + x * (-16.0)))) / 24.0;
    case 4:
      return (1.0 + x * (6.0 + x * (12.0 + x * (8.0)))) / 48.0;
    }
  case 5:
    switch (i) {
    case 0:
      return (-1.0 + x * (8.0 + x * (-24.0 + x * (32.0 + x * (-16))))) / 384.0;
    case 1:
      return (-75.0 + x * (168.0 + x * (-72.0 + x * (-96.0 + x * (80.0))))) /
             384.0;
    case 2:
      return (-77.0 + x * (-88.0 + x * (168.0 + x * (32.0 + x * (-80.0))))) /
             192.0;
    case 3:
      return (77.0 + x * (-88.0 + x * (-168.0 + x * (32.0 + x * (80.0))))) /
             192.0;
    case 4:
      return (75.0 + x * (168.0 + x * (72.0 + x * (-96.0 + x * (-80))))) /
             384.0;
    case 5:
      return (1.0 + x * (8.0 + x * (24.0 + x * (32.0 + x * (16.0))))) / 384.0;
    }
  case 6:
    switch (i) {
    case 0:
      return (-1.0 +
              x * (10.0 + x * (-40.0 + x * (80.0 + x * (-80.0 + x * 32.0))))) /
             3840.0;
    case 1:
      return (-59.0 +
              x * (185.0 + x * (-200.0 + x * (40.0 + x * (80.0 - x * 48.0))))) /
             960.0;
    case 2:
      return (-289.0 +
              x * (158.0 +
                   x * (344.0 + x * (-272.0 + x * (-80.0 + x * 96.0))))) /
             768.0;
    case 3:
      return (x * (-77.0 + x * x * (56.0 - x * x * 16.0))) / 96.0;
    case 4:
      return (289.0 +
              x * (158.0 +
                   x * (-344.0 + x * (-272.0 + x * (80.0 + x * 96.0))))) /
             768.0;
    case 5:
      return (59.0 +
              x * (185.0 + x * (200.0 + x * (40.0 + x * (-80.0 - x * 48.0))))) /
             960.0;
    case 6:
      return (1.0 +
              x * (10.0 + x * (40.0 + x * (80.0 + x * (80.0 + x * 32.0))))) /
             3840.0;
    }
  }

  throw std::runtime_error("Internal interpolation error.");
}
} // namespace Utils

#endif
