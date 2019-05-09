/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
     Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file
 *  Implementation of dpd.hpp.
 */
#include "dpd.hpp"

#ifdef DPD

#include "communication.hpp"
#include "random.hpp"
#include "short_range_loop.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
using Utils::Vector3d;
#include <utils/NoOp.hpp>
#include <utils/constants.hpp>
#include <utils/math/tensor_product.hpp>

#include <boost/mpi/collectives/reduce.hpp>

void dpd_heat_up() {
  double pref_scale = sqrt(3);
  dpd_update_params(pref_scale);
}

void dpd_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  dpd_update_params(pref_scale);
}

int dpd_set_params(int part_type_a, int part_type_b, double gamma, double r_c,
                   int wf, double tgamma, double tr_c, int twf) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  data->dpd_gamma = gamma;
  data->dpd_r_cut = r_c;
  data->dpd_wf = wf;
  data->dpd_pref2 = sqrt(24.0 * temperature * gamma / time_step);
  data->dpd_tgamma = tgamma;
  data->dpd_tr_cut = tr_c;
  data->dpd_twf = twf;
  data->dpd_pref4 = sqrt(24.0 * temperature * tgamma / time_step);

  /* Only make active if the DPD thermostat is
     activated, otherwise it will by activated
     by dpd_init() on thermostat change.
  */
  if (thermo_switch & THERMO_DPD) {
    data->dpd_pref1 = gamma / time_step;
    data->dpd_pref3 = tgamma / time_step;
  } else {
    data->dpd_pref1 = 0.0;
    data->dpd_pref3 = 0.0;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

void dpd_init() {
  for (int type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (int type_b = 0; type_b < max_seen_particle_type; type_b++) {
      auto data = get_ia_param(type_a, type_b);
      if ((data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0)) {
        data->dpd_pref1 = data->dpd_gamma / time_step;
        data->dpd_pref2 =
            sqrt(24.0 * temperature * data->dpd_gamma / time_step);
        data->dpd_pref3 = data->dpd_tgamma / time_step;
        data->dpd_pref4 =
            sqrt(24.0 * temperature * data->dpd_tgamma / time_step);
      }
    }
  }
}

void dpd_switch_off() {
  for (int type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (int type_b = 0; type_b < max_seen_particle_type; type_b++) {
      auto data = get_ia_param(type_a, type_b);
      data->dpd_pref1 = data->dpd_pref3 = 0.0;
    }
  }
}

void dpd_update_params(double pref_scale) {
  int type_a, type_b;
  IA_parameters *data;

  for (type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (type_b = 0; type_b < max_seen_particle_type; type_b++) {
      data = get_ia_param(type_a, type_b);
      if ((data->dpd_r_cut != 0) || (data->dpd_tr_cut != 0)) {
        data->dpd_pref2 *= pref_scale;
        data->dpd_pref4 *= pref_scale;
      }
    }
  }
}

static double weight(int type, double r_cut, double dist_inv) {
  if (type == 0) {
    return dist_inv;
  }
  return dist_inv - 1.0 / r_cut;
}

using Utils::Vector;

template <class T, size_t N, size_t M> using Matrix = Vector<Vector<T, M>, N>;

/**
 * @brief Kronecker delta.
 *
 * @return true iff I = J.
 *
 */
template <size_t I, size_t J> constexpr bool delta() { return I == J; }

namespace detail {
template <class T, size_t I, size_t... J>
constexpr Vector<T, sizeof...(J)> e_impl(std::index_sequence<J...>) {
  return {delta<I, J>()...};
}
} // namespace detail

/**
 * @brief I-th vector of the N-dimensional standard basis.
 * @tparam T Floating point type
 * @tparam N Dimension
 * @tparam I Index
 *                   .-> I-th element
 * @return { 0, ..., 1, 0, ... }
 */
template <class T, size_t N, size_t I> constexpr Vector<T, N> e() {
  return detail::e_impl<T, I>(std::make_index_sequence<N>{});
}

namespace detail {
template <class T, size_t... I>
constexpr Matrix<T, sizeof...(I), sizeof...(I)>
I_impl(std::index_sequence<I...>) {
  return {e<T, sizeof...(I), I>()...};
}
} // namespace detail

/**
 * @briefn NxN identity matrix
 * @tparam T Floating point type
 * @tparam N Dimension
 * @return e<T, I>() * e<T, J>()
 */
template <class T, size_t N> constexpr Matrix<T, N, N> I() {
  return detail::I_impl<T>(std::make_index_sequence<N>{});
}

Vector3d dpd_pair_force(Particle const *p1, Particle const *p2,
                        const IA_parameters *ia_params, double const *d,
                        double dist, double dist2) {
  using Utils::tensor_product;
  Vector3d f{};

  auto const dist_inv = 1.0 / dist;
  auto const v21 = p1->m.v - p2->m.v;
  auto const d21 = Utils::Vector3d{d[0], d[1], d[2]};

  /* Projection to d21 subspace */
  auto const P = Matrix<double, 3, 1>{{d[0], d[1], d[2]}};

  if ((dist < ia_params->dpd_r_cut) && (ia_params->dpd_gamma > 0.0)) {
    auto const omega =
        weight(ia_params->dpd_wf, ia_params->dpd_r_cut, dist_inv);
    auto const omega2 = Utils::sqr(omega);

    auto const friction = ia_params->dpd_gamma * omega2 * (P * v21)[0];
    double noise = (ia_params->dpd_pref2 > 0.0) ? (d_random() - 0.5) : 0.0;

    f += (ia_params->dpd_pref2 * omega * noise - friction) * d21;
  }
  // DPD2 part
  if ((dist < ia_params->dpd_tr_cut) && (ia_params->dpd_pref3 > 0.0)) {
    auto const omega =
        weight(ia_params->dpd_twf, ia_params->dpd_tr_cut, dist_inv);
    auto const omega2 = Utils::sqr(omega);

    auto const P = dist2 * I<double, 3>() - tensor_product(d21, d21);

    auto const noise_vec =
        (ia_params->dpd_pref2 > 0.0)
            ? Vector3d{d_random() - 0.5, d_random() - 0.5, d_random() - 0.5}
            : Vector3d{};

    auto f_D = P * v21;
    auto f_R = P * noise_vec;

    f_D *= ia_params->dpd_pref3 * omega2 * time_step;
    f_R *= ia_params->dpd_pref4 * omega * dist_inv;

    f += f_R - f_D;
  }

  return f;
}

static auto dpd_stress_local() {
  Utils::Vector<Utils::Vector3d, 3> stress{};

  short_range_loop(
      Utils::NoOp{},
      [&stress](const Particle &p1, const Particle &p2, Distance const &d) {
        auto const f =
            dpd_pair_force(&p1, &p2, get_ia_param_safe(p1.p.type, p2.p.type),
                           d.vec21.data(), sqrt(d.dist2), d.dist2);

        stress += Utils::tensor_product(d.vec21, f);
      });

  return stress;
}

void mpi_dpd_stress_slave() {
  boost::mpi::reduce(comm_cart, dpd_stress_local(), std::plus<>(), 0);
}
REGISTER_CALLBACK(mpi_dpd_stress_slave)

Utils::Vector9d mpi_dpd_stress() {
  Utils::Vector<Utils::Vector3d, 3> stress{};

  mpi_call(mpi_dpd_stress_slave);
  boost::mpi::reduce(comm_cart, dpd_stress_local(), stress, std::plus<>(), 0);

  return {stress[0][0], stress[0][1], stress[0][2], stress[1][0], stress[1][1],
          stress[1][2], stress[2][0], stress[2][1], stress[2][2]};
}

Utils::Vector9d dpd_stress() {
  return mpi_dpd_stress() / (box_l[0] * box_l[1] * box_l[2]);
}
#endif
