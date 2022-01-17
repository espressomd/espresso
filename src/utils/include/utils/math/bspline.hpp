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
#include "utils/device_qualifier.hpp"

#include <stdexcept>
#include <type_traits>

namespace Utils {
/** @brief Formula of the B-spline. */
template <int order, typename T>
DEVICE_QUALIFIER auto bspline(int i, T x)
    -> std::enable_if_t<(order > 0) && (order <= 7), T> {
  DEVICE_ASSERT(i < order);
  DEVICE_ASSERT(x >= T(-0.5));
  DEVICE_ASSERT(x <= T(0.5));

  switch (order) {
  case 1:
    return T(1.);
  case 2:
    switch (i) {
    case 0:
      return T(0.5) - x;
    case 1:
      return T(0.5) + x;
    }
  case 3:
    switch (i) {
    case 0:
      return T(0.5) * sqr(T(0.5) - x);
    case 1:
      return T(0.75) - sqr(x);
    case 2:
      return T(0.5) * sqr(T(0.5) + x);
    }
  case 4:
    switch (i) {
    case 0:
      return (T(1.) + x * (T(-6.) + x * (T(12.) - x * T(8.)))) / T(48.);
    case 1:
      return (T(23.) + x * (T(-30.) + x * (T(-12.) + x * T(24.)))) / T(48.);
    case 2:
      return (T(23.) + x * (T(30.) + x * (T(-12.) - x * T(24.)))) / T(48.);
    case 3:
      return (T(1.) + x * (T(6.) + x * (T(12.) + x * T(8.)))) / T(48.);
    }
  case 5:
    switch (i) {
    case 0:
      return (T(1.) +
              x * (T(-8.) + x * (T(24.) + x * (T(-32.) + x * T(16.))))) /
             T(384.);
    case 1:
      return (T(19.) +
              x * (T(-44.) + x * (T(24.) + x * (T(16.) - x * T(16.))))) /
             T(96.);
    case 2:
      return (T(115.) + x * x * (T(-120.) + x * x * T(48.))) / T(192.);
    case 3:
      return (T(19.) +
              x * (T(44.) + x * (T(24.) + x * (T(-16.) - x * T(16.))))) /
             T(96.);
    case 4:
      return (T(1.) + x * (T(8.) + x * (T(24.) + x * (T(32.) + x * T(16.))))) /
             T(384.);
    }
  case 6:
    switch (i) {
    case 0:
      return (T(1.) +
              x * (T(-10.) +
                   x * (T(40.) + x * (T(-80.) + x * (T(80.) - x * T(32.)))))) /
             T(3840.);
    case 1:
      return (T(237.) +
              x * (T(-750.) +
                   x * (T(840.) +
                        x * (T(-240.) + x * (T(-240.) + x * T(160.)))))) /
             T(3840.);
    case 2:
      return (T(841.) +
              x * (T(-770.) +
                   x * (T(-440.) +
                        x * (T(560.) + x * (T(80.) - x * T(160.)))))) /
             T(1920.);
    case 3:
      return (T(841.) +
              x * (+T(770.) +
                   x * (T(-440.) +
                        x * (T(-560.) + x * (T(80.) + x * T(160.)))))) /
             T(1920.);
    case 4:
      return (T(237.) +
              x * (T(750.) +
                   x * (T(840.) +
                        x * (T(240.) + x * (T(-240.) - x * T(160.)))))) /
             T(3840.);
    case 5:
      return (T(1.) +
              x * (T(10.) +
                   x * (T(40.) + x * (T(80.) + x * (T(80.) + x * T(32.)))))) /
             T(3840.);
    }
  case 7:
    switch (i) {
    case 0:
      return (T(1.) +
              x * (T(-12.) +
                   x * (T(60.) +
                        x * (T(-160.) +
                             x * (T(240.) + x * (T(-192.) + x * T(64.))))))) /
             T(46080.);
    case 1:
      return (T(361.) +
              x * (T(-1416.) +
                   x * (T(2220.) +
                        x * (T(-1600.) +
                             x * (T(240.) + x * (T(384.) - x * T(192.))))))) /
             T(23040.);
    case 2:
      return (T(10543.) +
              x * (T(-17340.) +
                   x * (T(4740.) +
                        x * (T(6880.) + x * (T(-4080.) +
                                             x * (T(-960.) + x * T(960.))))))) /
             T(46080.);
    case 3:
      return (T(5887.) +
              x * x * (T(-4620.) + x * x * (T(1680.) - x * x * T(320.)))) /
             T(11520.);
    case 4:
      return (T(10543.) +
              x * (T(17340.) +
                   x * (T(4740.) +
                        x * (T(-6880.) +
                             x * (T(-4080.) + x * (T(960.) + x * T(960.))))))) /
             T(46080.);
    case 5:
      return (T(361.) +
              x * (T(1416.) +
                   x * (T(2220.) +
                        x * (T(1600.) +
                             x * (T(240.) + x * (T(-384.) - x * T(192.))))))) /
             T(23040.);
    case 6:
      return (T(1.) +
              x * (T(12.) +
                   x * (T(60.) +
                        x * (T(160.) +
                             x * (T(240.) + x * (T(192.) + x * T(64.))))))) /
             T(46080.);
    }
  }

  DEVICE_THROW(std::runtime_error("Internal interpolation error."));
  return T{};
}

/**
 * @brief Calculate B-splines.
 * @param i knot number, using 0-based indexing
 * @param x position in the range (-0.5, 0.5)
 * @param k order of the B-spline, using 1-based indexing, i.e. a
 * B-spline of order @p k is a polynomial of degree <tt>k-1</tt>
 */
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

  return T(0.);
}

/** @brief Derivative of the B-spline. */
template <int order, typename T = double>
DEVICE_QUALIFIER auto bspline_d(int i, T x)
    -> std::enable_if_t<(order > 0) && (order <= 7), T> {
  DEVICE_ASSERT(i < order);
  DEVICE_ASSERT(x >= T(-0.5));
  DEVICE_ASSERT(x <= T(0.5));

  switch (order - 1) {
  case 0:
    return 0.;
  case 1:
    switch (i) {
    case 0:
      return T(-1.);
    case 1:
      return T(1.);
    }
  case 2:
    switch (i) {
    case 0:
      return x - T(0.5);
    case 1:
      return T(-2.) * x;
    case 2:
      return x + T(0.5);
    }
  case 3:
    switch (i) {
    case 0:
      return (T(-1.) + x * (T(4.) + x * T(-4.))) / T(8.);
    case 1:
      return (T(-5.) + x * (T(-4.) + x * T(12.))) / T(8.);
    case 2:
      return (T(5.) + x * (T(-4.) + x * T(-12.))) / T(8.);
    case 3:
      return (T(1.) + x * (T(4.) + x * T(4.))) / T(8.);
    }
  case 4:
    switch (i) {
    case 0:
      return (T(-1.) + x * (T(6.) + x * (T(-12.) + x * T(8.)))) / T(48.);
    case 1:
      return (T(-11.) + x * (T(12.) + x * (T(12.) + x * T(-16.)))) / T(24.);
    case 2:
      return (x * (T(-5.) + x * x * T(4.))) / T(4.);
    case 3:
      return (T(11.) + x * (T(12.) + x * (T(-12.) + x * T(-16.)))) / T(24.);
    case 4:
      return (T(1.) + x * (T(6.) + x * (T(12.) + x * T(8.)))) / T(48.);
    }
  case 5:
    switch (i) {
    case 0:
      return (T(-1.) +
              x * (T(8.) + x * (T(-24.) + x * (T(32.) + x * T(-16.))))) /
             T(384.);
    case 1:
      return (T(-75.) +
              x * (T(168.) + x * (T(-72.) + x * (T(-96.) + x * (T(80.)))))) /
             T(384.);
    case 2:
      return (T(-77.) +
              x * (T(-88.) + x * (T(168.) + x * (T(32.) + x * T(-80.))))) /
             T(192.);
    case 3:
      return (T(77.) +
              x * (T(-88.) + x * (T(-168.) + x * (T(32.) + x * T(80.))))) /
             T(192.);
    case 4:
      return (T(75.) +
              x * (T(168.) + x * (T(72.) + x * (T(-96.) + x * T(-80.))))) /
             T(384.);
    case 5:
      return (T(1.) + x * (T(8.) + x * (T(24.) + x * (T(32.) + x * T(16.))))) /
             T(384.);
    }
  case 6:
    switch (i) {
    case 0:
      return (T(-1.) +
              x * (T(10.) +
                   x * (T(-40.) + x * (T(80.) + x * (T(-80.) + x * T(32.)))))) /
             T(3840.);
    case 1:
      return (T(-59.) +
              x * (T(185.) +
                   x * (T(-200.) + x * (T(40.) + x * (T(80.) - x * T(48.)))))) /
             T(960.);
    case 2:
      return (T(-289.) +
              x * (T(158.) +
                   x * (T(344.) +
                        x * (T(-272.) + x * (T(-80.) + x * T(96.)))))) /
             T(768.);
    case 3:
      return (x * (T(-77.) + x * x * (T(56.) - x * x * T(16.)))) / T(96.);
    case 4:
      return (T(289.) + x * (T(158.) + x * (T(-344.) +
                                            x * (T(-272.) +
                                                 x * (T(80.) + x * T(96.)))))) /
             T(768.);
    case 5:
      return (T(59.) +
              x * (T(185.) +
                   x * (T(200.) + x * (T(40.) + x * (T(-80.) - x * T(48.)))))) /
             T(960.);
    case 6:
      return (T(1.) +
              x * (T(10.) +
                   x * (T(40.) + x * (T(80.) + x * (T(80.) + x * T(32.)))))) /
             T(3840.);
    }
  }

  DEVICE_THROW(std::runtime_error("Internal interpolation error."));
  return T{};
}
} // namespace Utils

#endif
