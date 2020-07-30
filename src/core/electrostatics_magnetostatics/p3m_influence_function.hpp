/*
 * Copyright (C) 2019-2020 The ESPResSo project
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
#ifndef ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP
#define ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/int_pow.hpp>
#include <utils/math/sinc.hpp>
#include <utils/math/sqr.hpp>

#include <boost/range/numeric.hpp>

#include <cmath>

/**
 * @brief Functor for calculating the Hockney/Eastwood/Ballenegger
 *        generic P3M influence function.
 * @tparam cao Charge assignment order.
 * @tparam S Power of the differential operator e.g. 0 for energy,
 *           1 for force and so on.
 * @tparam m Number of aliasing terms to take into account.
 */
template <size_t cao, size_t S, size_t m> struct InfluenceFunction {
private:
  enum : int { RX = 0, RY = 1, RZ = 2 };
  enum : int { KY = 0, KZ = 1, KX = 2 };

  std::array<std::vector<int>, 3>
  calc_meshift(std::array<int, 3> const &mesh_size) const {
    std::array<std::vector<int>, 3> ret;

    for (size_t i = 0; i < 3; i++) {
      ret[i].resize(mesh_size[i]);

      ret[i][0] = 0;
      for (int j = 1; j <= mesh_size[i] / 2; j++) {
        ret[i][j] = j;
        ret[i][mesh_size[i] - j] = -j;
      }
    }

    return ret;
  }

  template <typename T> T g_ewald(T alpha, T k2) const {
    auto constexpr limit = T{30};
    auto const exponent = Utils::sqr(1. / (2. * alpha)) * k2;
    return (exponent < limit) ? std::exp(-exponent) * (4. * Utils::pi<T>() / k2)
                              : 0.0;
  }

  std::pair<double, double> aliasing_sums_ik(double alpha,
                                             const Utils::Vector3d &k,
                                             const Utils::Vector3d &h) const {
    using Utils::int_pow;
    using Utils::pi;
    using Utils::sinc;
    using Utils::Vector3d;

    double numerator = 0.0;
    double denominator = 0.0;

    for (int mx = -m; mx <= m; mx++) {
      for (int my = -m; my <= m; my++) {
        for (int mz = -m; mz <= m; mz++) {
          auto const km =
              k + 2 * pi() * Vector3d{mx / h[RX], my / h[RY], mz / h[RZ]};
          auto const U2 = int_pow<2 * cao>(sinc(0.5 * km[RX] * h[RX] / pi()) *
                                           sinc(0.5 * km[RY] * h[RY] / pi()) *
                                           sinc(0.5 * km[RZ] * h[RZ] / pi()));

          numerator += U2 * g_ewald(alpha, km.norm2()) * int_pow<S>(k * km);
          denominator += U2;
        }
      }
    }
    return {numerator, denominator};
  }

  double G_opt(double alpha, const Utils::Vector3d &k,
               const Utils::Vector3d &h) const {
    using Utils::int_pow;
    using Utils::sqr;

    auto const k2 = k.norm2();
    if (k2 == 0.0) {
      return 0.0;
    }

    const auto as = aliasing_sums_ik(alpha, k, h);
    return as.first / (int_pow<S>(k2) * sqr(as.second));
  }

public:
  std::vector<double> operator()(const P3MParameters &params,
                                 const fft_data_struct &fft,
                                 const Utils::Vector3d &box_l) const {
    auto const shifts =
        calc_meshift({params.mesh[0], params.mesh[1], params.mesh[2]});

    auto const size =
        boost::accumulate(fft.plan[3].new_mesh, 1, std::multiplies<>());
    auto const start = Utils::Vector3i{fft.plan[3].start};
    auto const end = start + Utils::Vector3i{fft.plan[3].new_mesh};

    /* The influence function grid */
    auto g = std::vector<double>(size, 0.);

    /* Skip influence function calculation in tuning mode,
       the results need not be correct for timing. */
    if (params.tuning) {
      return g;
    }

    auto const h = Utils::Vector3d{params.a};

    Utils::Vector3i n{};
    for (n[0] = fft.plan[3].start[0]; n[0] < end[0]; n[0]++) {
      for (n[1] = fft.plan[3].start[1]; n[1] < end[1]; n[1]++) {
        for (n[2] = fft.plan[3].start[2]; n[2] < end[2]; n[2]++) {
          auto const ind = Utils::get_linear_index(
              n - start, Utils::Vector3i{fft.plan[3].new_mesh},
              Utils::MemoryOrder::ROW_MAJOR);

          if ((n[KX] % (params.mesh[RX] / 2) == 0) &&
              (n[KY] % (params.mesh[RY] / 2) == 0) &&
              (n[KZ] % (params.mesh[RZ] / 2) == 0)) {
            g[ind] = 0.0;
          } else {
            auto const k = 2 * Utils::pi() *
                           Utils::Vector3d{shifts[RX][n[KX]] / box_l[RX],
                                           shifts[RY][n[KY]] / box_l[RY],
                                           shifts[RZ][n[KZ]] / box_l[RZ]};

            g[ind] = G_opt(params.alpha, k, h);
          }
        }
      }
    }

    return g;
  }
};

#endif // ESPRESSO_P3M_INFLUENCE_FUNCTION_HPP
