/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "electrostatics/elc.hpp"

#include "electrostatics/coulomb.hpp"
#include "electrostatics/elc_parallel.hpp"
#include "electrostatics/p3m.hpp"

#include "BoxGeometry.hpp"
#include "aosoa_pack.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "kokkos_helpers.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <Kokkos_Core.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <omp.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

/** ELC charge sum/assign protocol: real charges, image charges, or both. */
enum class ChargeProtocol : int { REAL, IMAGE, BOTH };

namespace {

enum class FreqKind : int { P, Q, PQ };

/** @brief One term of the far-field sum. */
struct FrequencyEntry {
  FreqKind kind;
  std::size_t p; ///< x frequency index (>= 1 for P and PQ)
  std::size_t q; ///< y frequency index (>= 1 for Q and PQ)
  double omega;
  double pref;    ///< scales the real and layer-image block
  double pref_di; ///< scales the image-sum block (dielectric contrast only)
  double fac_delta_mid_bot;
  double fac_delta_mid_top;
  double fac_delta;
  std::size_t k;            ///< 4 for P and Q, 8 for PQ
  std::size_t sums_offset;  ///< into the pass-1 buffer (k or 2k slots)
  std::size_t coeff_offset; ///< into the coefficient buffer (k slots)
};

struct FarFieldTable {
  std::vector<FrequencyEntry> entries; ///< all P, then all Q, then all PQ
  std::array<std::size_t, 3> n_kind{}; ///< number of P, Q and PQ entries
  std::size_t n_sums = 0;              ///< total pass-1 slots
  std::size_t n_coeffs = 0;            ///< total coefficient slots
};

template <FreqKind kind>
using FreqKindTag = std::integral_constant<FreqKind, kind>;

/**
 * @brief Call @p kernel on every entry in table order, with the frequency
 * kind as a compile-time tag, so that the kernel body has fixed sizes.
 */
template <class Kernel>
inline void for_each_entry(FarFieldTable const &table, Kernel &&kernel) {
  auto const *entry = table.entries.data();
  for (std::size_t n = 0; n < table.n_kind[0]; ++n) {
    kernel(*entry++, FreqKindTag<FreqKind::P>{});
  }
  for (std::size_t n = 0; n < table.n_kind[1]; ++n) {
    kernel(*entry++, FreqKindTag<FreqKind::Q>{});
  }
  for (std::size_t n = 0; n < table.n_kind[2]; ++n) {
    kernel(*entry++, FreqKindTag<FreqKind::PQ>{});
  }
}

std::pair<std::size_t, std::size_t> n_frequencies(BoxGeometry const &box_geo,
                                                  double far_cut) {
  assert(far_cut >= 0.);
  auto const n_freq_x =
      static_cast<std::size_t>(std::ceil(far_cut * box_geo.length()[0]) + 1.);
  auto const n_freq_y =
      static_cast<std::size_t>(std::ceil(far_cut * box_geo.length()[1]) + 1.);
  return {n_freq_x, n_freq_y};
}

FrequencyEntry make_entry(elc_data const &elc, double prefactor,
                          BoxGeometry const &box_geo, FreqKind kind,
                          std::size_t p, std::size_t q, double omega) {
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const k = (kind == FreqKind::PQ) ? std::size_t{8u} : std::size_t{4u};
  auto const pref_di = prefactor * ((kind == FreqKind::PQ) ? 8. : 4.) *
                       std::numbers::pi * xy_area_inv;
  auto const pref = -pref_di / expm1(omega * box_geo.length()[2]);
  FrequencyEntry entry{kind, p, q, omega, pref, pref_di, 1., 1., 1., k, 0, 0};
  if (elc.dielectric_contrast_on) {
    auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
    auto const fac_elc = 1. / (1. - delta * exp(-omega * 2. * elc.box_h));
    entry.fac_delta_mid_bot = elc.delta_mid_bot * fac_elc;
    entry.fac_delta_mid_top = elc.delta_mid_top * fac_elc;
    entry.fac_delta = entry.fac_delta_mid_bot * elc.delta_mid_top;
  }
  return entry;
}

/** @brief List the far-field terms in the order of the original loops. */
FarFieldTable make_frequency_table(elc_data const &elc, double prefactor,
                                   BoxGeometry const &box_geo) {
  auto constexpr c_2pi = 2. * std::numbers::pi;
  auto const [n_freq_x, n_freq_y] = n_frequencies(box_geo, elc.far_cut);
  auto const ux = box_geo.length_inv()[0];
  auto const uy = box_geo.length_inv()[1];
  FarFieldTable table;
  /* the second condition is just for the case of numerical accident */
  for (std::size_t p = 1;
       ux * static_cast<double>(p - 1) < elc.far_cut && p <= n_freq_x; p++) {
    auto const omega = c_2pi * ux * static_cast<double>(p);
    table.entries.push_back(
        make_entry(elc, prefactor, box_geo, FreqKind::P, p, 0, omega));
  }
  for (std::size_t q = 1;
       uy * static_cast<double>(q - 1) < elc.far_cut && q <= n_freq_y; q++) {
    auto const omega = c_2pi * uy * static_cast<double>(q);
    table.entries.push_back(
        make_entry(elc, prefactor, box_geo, FreqKind::Q, 0, q, omega));
  }
  for (std::size_t p = 1;
       ux * static_cast<double>(p - 1) < elc.far_cut && p <= n_freq_x; p++) {
    for (std::size_t q = 1;
         Utils::sqr(ux * static_cast<double>(p - 1)) +
                 Utils::sqr(uy * static_cast<double>(q - 1)) <
             elc.far_cut2 &&
         q <= n_freq_y;
         q++) {
      auto const omega = c_2pi * sqrt(Utils::sqr(ux * static_cast<double>(p)) +
                                      Utils::sqr(uy * static_cast<double>(q)));
      table.entries.push_back(
          make_entry(elc, prefactor, box_geo, FreqKind::PQ, p, q, omega));
    }
  }
  auto const blocks = elc.dielectric_contrast_on ? 2u : 1u;
  for (auto &entry : table.entries) {
    ++table.n_kind[static_cast<std::size_t>(entry.kind)];
    entry.sums_offset = table.n_sums;
    entry.coeff_offset = table.n_coeffs;
    table.n_sums += blocks * entry.k;
    table.n_coeffs += entry.k;
  }
  return table;
}

using HostMatrix =
    Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace>;

/** @brief Sin/cos of @f$ 2\pi f u x_{dir} @f$, one row per particle. */
struct SinCosCache {
  HostMatrix sin;
  HostMatrix cos;
};

SinCosCache make_sc_cache(auto const &aosoa, std::size_t n_part,
                          std::size_t dir, std::size_t n_freq, double u) {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  auto constexpr c_2pi = 2. * std::numbers::pi;
  SinCosCache cache{
      HostMatrix(Kokkos::view_alloc("ELC::sin", Kokkos::WithoutInitializing),
                 n_part, n_freq),
      HostMatrix(Kokkos::view_alloc("ELC::cos", Kokkos::WithoutInitializing),
                 n_part, n_freq)};
  kokkos_parallel_range_for<execution_space>(
      "ELC::sc_cache", std::size_t{0u}, n_part,
      [sin_view = cache.sin, cos_view = cache.cos, &aosoa, dir, n_freq,
       u](std::size_t i) {
        auto const x = aosoa.position(i, dir);
        for (std::size_t freq = 1; freq <= n_freq; freq++) {
          auto const arg = c_2pi * u * static_cast<double>(freq) * x;
          sin_view(i, freq - 1) = sin(arg);
          cos_view(i, freq - 1) = cos(arg);
        }
      });
  return cache;
}

/** @brief Number of trigonometric products per particle and frequency. */
template <FreqKind kind>
inline constexpr std::size_t n_trig = (kind == FreqKind::PQ) ? 4u : 2u;

template <FreqKind kind> using TrigTerms = std::array<double, n_trig<kind>>;
template <FreqKind kind>
using RealBlock = std::array<double, 2u * n_trig<kind>>;

/** @brief s, c (P, Q) or ss, sc, cs, cc (PQ) of one particle. */
template <FreqKind kind>
inline TrigTerms<kind> particle_terms(FrequencyEntry const &entry,
                                      SinCosCache const &scx,
                                      SinCosCache const &scy, std::size_t i) {
  if constexpr (kind == FreqKind::P) {
    return {scx.sin(i, entry.p - 1), scx.cos(i, entry.p - 1)};
  } else if constexpr (kind == FreqKind::Q) {
    return {scy.sin(i, entry.q - 1), scy.cos(i, entry.q - 1)};
  } else {
    auto const sx = scx.sin(i, entry.p - 1);
    auto const cx = scx.cos(i, entry.p - 1);
    auto const sy = scy.sin(i, entry.q - 1);
    auto const cy = scy.cos(i, entry.q - 1);
    return {sx * sy, sx * cy, cx * sy, cx * cy};
  }
}

/** @brief Real-charge block of one particle (P half, then M half). */
template <FreqKind kind>
inline RealBlock<kind> real_block(TrigTerms<kind> const &t, double q,
                                  double e) {
  auto constexpr h = n_trig<kind>;
  RealBlock<kind> out;
  for (std::size_t m = 0; m < h; ++m) {
    if constexpr (kind == FreqKind::PQ) {
      out[m] = t[m] * q * e;
      out[h + m] = t[m] * q / e;
    } else {
      out[m] = q * t[m] * e;
      out[h + m] = q * t[m] / e;
    }
  }
  return out;
}

/**
 * @brief Pass 1: accumulate the sums of every frequency, reduce them over
 * threads and MPI ranks, and apply the prefactors.
 * @return @c table.n_coeffs coefficients, @c entry.k per entry
 */
std::vector<double>
far_field_coefficients(elc_data const &elc, FarFieldTable const &table,
                       SinCosCache const &scx, SinCosCache const &scy,
                       auto const &aosoa, std::size_t n_part) {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  ELC::PerThreadRows rows(table.n_sums);
  auto const contrast = elc.dielectric_contrast_on;
  auto const h = elc.box_h;
  kokkos_parallel_range_for<execution_space>(
      "ELC::far_field_sums", std::size_t{0u}, n_part, [&](std::size_t i) {
        auto *row = rows.row(omp_get_thread_num());
        auto const z = aosoa.position(i, 2);
        auto const q = aosoa.charge(i);
        for_each_entry(table, [&](FrequencyEntry const &entry, auto tag) {
          auto constexpr kind = decltype(tag)::value;
          auto constexpr half = n_trig<kind>;
          auto const omega = entry.omega;
          auto const t = particle_terms<kind>(entry, scx, scy, i);
          auto *real = row + entry.sums_offset;
          auto const blk = real_block<kind>(t, q, exp(omega * z));
          for (std::size_t m = 0; m < 2u * half; ++m) {
            real[m] += blk[m];
          }
          if (not contrast) {
            return;
          }
          auto *image = real + 2u * half;
          double e_p, e_m;
          if (z < elc.space_layer) { // handle the lower case first
            // negative sign is okay here as the image is located at -z
            auto const e = exp(-omega * z);
            auto const scale = q * elc.delta_mid_bot;
            for (std::size_t m = 0; m < half; ++m) {
              real[m] += scale * (t[m] * e);
              real[half + m] += scale * (t[m] / e);
            }
            e_p = (exp(omega * (-z - 2. * h)) * elc.delta_mid_bot +
                   exp(omega * (+z - 2. * h))) *
                  entry.fac_delta;
          } else {
            e_p = (exp(-omega * z) +
                   exp(omega * (z - 2. * h)) * elc.delta_mid_top) *
                  entry.fac_delta_mid_bot;
          }
          if (z > (h - elc.space_layer)) { // handle the upper case now
            auto const e = exp(omega * (2. * h - z));
            auto const scale = q * elc.delta_mid_top;
            for (std::size_t m = 0; m < half; ++m) {
              real[m] += scale * (t[m] * e);
              real[half + m] += scale * (t[m] / e);
            }
            e_m = (exp(omega * (+z - 4. * h)) * elc.delta_mid_top +
                   exp(omega * (-z - 2. * h))) *
                  entry.fac_delta;
          } else {
            e_m = (exp(omega * (+z - 2. * h)) +
                   exp(omega * (-z - 2. * h)) * elc.delta_mid_bot) *
                  entry.fac_delta_mid_top;
          }
          for (std::size_t m = 0; m < half; ++m) {
            if constexpr (kind == FreqKind::PQ) {
              image[m] += t[m] * (e_p * q);
              image[half + m] += t[m] * (e_m * q);
            } else {
              image[m] += q * t[m] * e_p;
              image[half + m] += q * t[m] * e_m;
            }
          }
        });
      });
  Kokkos::fence();

  auto const local = rows.fold();
  std::vector<double> global(table.n_sums);
  boost::mpi::all_reduce(comm_cart, local.data(),
                         static_cast<int>(table.n_sums), global.data(),
                         std::plus<>());

  std::vector<double> coeffs(table.n_coeffs);
  for (auto const &entry : table.entries) {
    for (std::size_t m = 0; m < entry.k; ++m) {
      auto value = entry.pref * global[entry.sums_offset + m];
      if (contrast) {
        value += entry.pref_di * global[entry.sums_offset + entry.k + m];
      }
      coeffs[entry.coeff_offset + m] = value;
    }
  }
  return coeffs;
}
} // namespace

void ElectrostaticLayerCorrection::check_gap(CellStructure const &cs) const {
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  auto const box_h = elc.box_h;
  auto const in_gap = [&aosoa, box_h](std::size_t i) {
    auto const z = aosoa.position(i, 2);
    return aosoa.charge(i) != 0. and (z < 0. or z > box_h);
  };
  auto const n_in_gap = ELC::parallel_array_sum<1>(
      "ELC::check_gap", n_part,
      [&in_gap](std::size_t i, std::array<double, 1> &acc) {
        if (in_gap(i)) {
          acc[0] += 1.;
        }
      })[0];
  if (n_in_gap == 0.) {
    return;
  }
  // runtimeErrorMsg() is not thread-safe: report serially
  for (std::size_t i = 0; i < n_part; ++i) {
    if (in_gap(i)) {
      auto const z = aosoa.position(i, 2);
      runtimeErrorMsg() << "Particle " << aosoa.id(i) << " entered ELC gap "
                        << "region by " << ((z < 0.) ? z : z - box_h);
    }
  }
}

/*****************************************************************/
/* dipole terms */
/*****************************************************************/

/** Calculate the dipole force.
 *  See @cite yeh99a.
 */
void ElectrostaticLayerCorrection::add_dipole_force() const {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  auto const pref = prefactor * 4. * std::numbers::pi / box_geo.volume();

  /* for non-neutral systems, this shift gives the background contribution
   * (rsp. for this shift, the DM of the background is zero) */
  auto const shift = box_geo.length_half()[2];
  auto const mid = 0.5 * elc.box_h;

  check_gap(cs);

  // collect moments:
  // [0] sum q_i (z_i - L/2), [1] sum q_i (z_i - box_h/2), [2] sum q_i
  auto local = ELC::parallel_array_sum<3>(
      "ELC::dipole_force_moments", n_part,
      [this, &aosoa, shift, mid](std::size_t i, std::array<double, 3> &acc) {
        auto const q = aosoa.charge(i);
        auto const z = aosoa.position(i, 2);
        acc[0] += q * (z - shift);
        acc[1] += q * (z - mid);
        acc[2] += q;
        if (elc.dielectric_contrast_on) {
          if (z < elc.space_layer) {
            acc[0] += elc.delta_mid_bot * q * (-z - shift);
            acc[2] += elc.delta_mid_bot * q;
          }
          if (z > (elc.box_h - elc.space_layer)) {
            acc[0] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
            acc[2] += elc.delta_mid_top * q;
          }
        }
      });

  local[0] *= pref;
  local[1] *= pref / elc.box_h * box_geo.length()[2];
  local[2] *= pref;

  std::array<double, 3> global{};
  boost::mpi::all_reduce(comm_cart, local.data(), 3, global.data(),
                         std::plus<>());

  // Yeh + Berkowitz dipole term @cite yeh99a
  auto field_tot = global[0];

  // Constant potential contribution
  if (elc.const_pot) {
    auto const field_induced = global[1];
    auto const field_applied = elc.pot_diff / elc.box_h;
    field_tot -= field_applied + field_induced;
  }

  auto const neutralize = elc.neutralize;
  auto const background = global[2];
  auto scatter_force = cs.get_scatter_force();
  kokkos_parallel_range_for<execution_space>(
      "ELC::dipole_force", std::size_t{0u}, n_part,
      [&aosoa, &scatter_force, field_tot, neutralize, background,
       shift](std::size_t i) {
        auto const q = aosoa.charge(i);
        auto fz = -field_tot * q;
        if (!neutralize) {
          // SUBTRACT the forces of the P3M homogeneous neutralizing background
          fz += background * q * (aosoa.position(i, 2) - shift);
        }
        auto access = scatter_force.access();
        access(i, 2) += fz;
      });
}

/** Calculate the dipole energy.
 *  See @cite yeh99a.
 */
double ElectrostaticLayerCorrection::dipole_energy() const {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const pref = prefactor * 2. * std::numbers::pi / box_geo.volume();
  auto const lz = box_geo.length()[2];
  /* for nonneutral systems, this shift gives the background contribution
     (rsp. for this shift, the DM of the background is zero) */
  auto const shift = box_geo.length_half()[2];

  auto const &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();

  check_gap(cs);

  // [0] sum q_i               primary box
  // [1] sum q_i               boundary layers
  // [2] sum q_i (z_i - L/2)   primary box
  // [3] sum q_i (z_i - L/2)   boundary layers
  // [4] sum q_i (z_i - L/2)^2 primary box
  // [5] sum q_i (z_i - L/2)^2 boundary layers
  // [6] sum q_i z_i           primary box
  auto const local = ELC::parallel_array_sum<7>(
      "ELC::dipole_energy_moments", n_part,
      [this, &aosoa, shift](std::size_t i, std::array<double, 7> &acc) {
        auto const q = aosoa.charge(i);
        auto const z = aosoa.position(i, 2);
        acc[0] += q;
        acc[2] += q * (z - shift);
        acc[4] += q * (Utils::sqr(z - shift));
        acc[6] += q * z;
        if (elc.dielectric_contrast_on) {
          if (z < elc.space_layer) {
            acc[1] += elc.delta_mid_bot * q;
            acc[3] += elc.delta_mid_bot * q * (-z - shift);
            acc[5] += elc.delta_mid_bot * q * (Utils::sqr(-z - shift));
          }
          if (z > (elc.box_h - elc.space_layer)) {
            acc[1] += elc.delta_mid_top * q;
            acc[3] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
            acc[5] += elc.delta_mid_top * q *
                      (Utils::sqr(2. * elc.box_h - z - shift));
          }
        }
      });
  std::array<double, 7> global{};
  boost::mpi::all_reduce(comm_cart, local.data(), 7, global.data(),
                         std::plus<>());

  // Yeh + Berkowitz term @cite yeh99a
  auto energy = 2. * pref * (Utils::sqr(global[2]) + global[2] * global[3]);

  if (!elc.neutralize) {
    // SUBTRACT the energy of the P3M homogeneous neutralizing background
    energy +=
        2. * pref *
        (-global[0] * global[4] - (.25 - .5 / 3.) * Utils::sqr(global[0] * lz));
  }

  if (elc.dielectric_contrast_on) {
    if (elc.const_pot) {
      // zero potential difference contribution
      energy -= 2. * pref / elc.box_h * lz * global[6] *
                (global[6] - elc.box_h * global[0]);
      // external potential shift contribution
      energy -= 2. * elc.pot_diff / elc.box_h * global[6];
    }

    /* counter the P3M homogeneous background contribution to the
       boundaries. We never need that, since a homogeneous background
       spanning the artificial boundary layers is aphysical. */
    energy += pref * (-(global[1] * global[4] + global[0] * global[5]) -
                      (.5 - 1. / 3.) * global[0] * global[1] * Utils::sqr(lz));
  }

  return this_node == 0 ? energy : 0.;
}

/*****************************************************************/

struct ImageSum {
  double delta;
  double shift;
  double h;   // plate separation
  double dci; // delta complement inverse

  ImageSum(double delta, double shift, double h)
      : delta{delta}, shift{shift}, h{h}, dci{1. / (1. - delta)} {}

  /** @brief Image sum from the bottom layer. */
  double b(double q, double z) const {
    return q * dci * (z - 2. * delta * h * dci) - q * dci * shift;
  }

  /** @brief Image sum from the top layer. */
  double t(double q, double z) const {
    return q * dci * (z + 2. * delta * h * dci) - q * dci * shift;
  }
};

double ElectrostaticLayerCorrection::z_energy() const {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref = prefactor * 2. * std::numbers::pi * xy_area_inv;

  /* for non-neutral systems, this shift gives the background contribution
   * (rsp. for this shift, the DM of the background is zero) */
  auto const shift = box_geo.length_half()[2];

  std::array<double, 4> local{};
  if (elc.dielectric_contrast_on) {
    if (elc.const_pot) {
      // metallic boundaries
      local = ELC::parallel_array_sum<4>(
          "ELC::z_energy_metallic", n_part,
          [this, &aosoa, shift](std::size_t i, std::array<double, 4> &acc) {
            auto const z = aosoa.position(i, 2);
            auto const q = aosoa.charge(i);
            acc[0] += q;
            acc[1] += q * (z - shift);
            if (z < elc.space_layer) {
              acc[2] -= elc.delta_mid_bot * q;
              acc[3] -= elc.delta_mid_bot * q * (-z - shift);
            }
            if (z > (elc.box_h - elc.space_layer)) {
              acc[2] += elc.delta_mid_top * q;
              acc[3] += elc.delta_mid_top * q * (2. * elc.box_h - z - shift);
            }
          });
    } else {
      // dielectric boundaries
      auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
      auto const fac_delta_mid_bot = elc.delta_mid_bot / (1. - delta);
      auto const fac_delta_mid_top = elc.delta_mid_top / (1. - delta);
      auto const fac_delta = delta / (1. - delta);
      auto const h = elc.box_h;
      ImageSum const image_sum{delta, shift, h};
      local = ELC::parallel_array_sum<4>(
          "ELC::z_energy_dielectric", n_part,
          [&](std::size_t i, std::array<double, 4> &acc) {
            auto const z = aosoa.position(i, 2);
            auto const q = aosoa.charge(i);
            acc[0] += q;
            acc[1] += q * (z - shift);
            if (z < elc.space_layer) {
              acc[2] += fac_delta * (elc.delta_mid_bot + 1.) * q;
              acc[3] +=
                  q * (image_sum.b(elc.delta_mid_bot * delta, -(2. * h + z)) +
                       image_sum.b(delta, -(2. * h - z)));
            } else {
              acc[2] += fac_delta_mid_bot * (1. + elc.delta_mid_top) * q;
              acc[3] += q * (image_sum.b(elc.delta_mid_bot, -z) +
                             image_sum.b(delta, -(2. * h - z)));
            }
            if (z > (h - elc.space_layer)) {
              // note the minus sign here which is required due to |z_i-z_j|
              acc[2] -= fac_delta * (elc.delta_mid_top + 1.) * q;
              acc[3] -=
                  q * (image_sum.t(elc.delta_mid_top * delta, 4. * h - z) +
                       image_sum.t(delta, 2. * h + z));
            } else {
              // note the minus sign here which is required due to |z_i-z_j|
              acc[2] -= fac_delta_mid_top * (1. + elc.delta_mid_bot) * q;
              acc[3] -= q * (image_sum.t(elc.delta_mid_top, 2. * h - z) +
                             image_sum.t(delta, 2. * h + z));
            }
          });
    }
  }
  std::array<double, 4> global{};
  boost::mpi::all_reduce(comm_cart, local.data(), 4, global.data(),
                         std::plus<>());

  auto const energy = global[1] * global[2] - global[0] * global[3];
  return (this_node == 0) ? -pref * energy : 0.;
}

void ElectrostaticLayerCorrection::add_z_force() const {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  if (not elc.dielectric_contrast_on) {
    return;
  }
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  auto const xy_area_inv = box_geo.length_inv()[0] * box_geo.length_inv()[1];
  auto const pref = prefactor * 2. * std::numbers::pi * xy_area_inv;

  std::array<double, 1> local{};
  if (elc.const_pot) {
    // metallic boundaries: just counter the 2 pi |z| contribution from P3M
    local = ELC::parallel_array_sum<1>(
        "ELC::z_force_metallic", n_part,
        [this, &aosoa](std::size_t i, std::array<double, 1> &acc) {
          auto const z = aosoa.position(i, 2);
          auto const q = aosoa.charge(i);
          if (z < elc.space_layer)
            acc[0] -= elc.delta_mid_bot * q;
          if (z > (elc.box_h - elc.space_layer))
            acc[0] += elc.delta_mid_top * q;
        });
  } else {
    // dielectric boundaries
    auto const delta = elc.delta_mid_top * elc.delta_mid_bot;
    auto const fac_delta_mid_bot = elc.delta_mid_bot / (1. - delta);
    auto const fac_delta_mid_top = elc.delta_mid_top / (1. - delta);
    auto const fac_delta = delta / (1. - delta);
    local = ELC::parallel_array_sum<1>(
        "ELC::z_force_dielectric", n_part,
        [&](std::size_t i, std::array<double, 1> &acc) {
          auto const z = aosoa.position(i, 2);
          auto const q = aosoa.charge(i);
          if (z < elc.space_layer) {
            acc[0] += fac_delta * (elc.delta_mid_bot + 1.) * q;
          } else {
            acc[0] += fac_delta_mid_bot * (elc.delta_mid_top + 1.) * q;
          }
          if (z > (elc.box_h - elc.space_layer)) {
            // note the minus sign here which is required due to |z_i-z_j|
            acc[0] -= fac_delta * (elc.delta_mid_top + 1.) * q;
          } else {
            // note the minus sign here which is required due to |z_i-z_j|
            acc[0] -= fac_delta_mid_top * (elc.delta_mid_bot + 1.) * q;
          }
        });
  }
  local[0] *= pref;
  double field = 0.;
  boost::mpi::all_reduce(comm_cart, local[0], field, std::plus<>());

  auto scatter_force = cs.get_scatter_force();
  kokkos_parallel_range_for<execution_space>(
      "ELC::z_force", std::size_t{0u}, n_part,
      [&aosoa, &scatter_force, field](std::size_t i) {
        auto access = scatter_force.access();
        access(i, 2) += field * aosoa.charge(i);
      });
}

/*****************************************************************/
/* far-field sums */
/*****************************************************************/

void ElectrostaticLayerCorrection::add_force() const {
  using execution_space = Kokkos::DefaultHostExecutionSpace;
  auto constexpr c_2pi = 2. * std::numbers::pi;
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();

  add_dipole_force();
  add_z_force();

  auto const table = make_frequency_table(elc, prefactor, box_geo);
  auto const [n_freq_x, n_freq_y] = n_frequencies(box_geo, elc.far_cut);
  auto const scx =
      make_sc_cache(aosoa, n_part, 0, n_freq_x, box_geo.length_inv()[0]);
  auto const scy =
      make_sc_cache(aosoa, n_part, 1, n_freq_y, box_geo.length_inv()[1]);
  auto const coeffs =
      far_field_coefficients(elc, table, scx, scy, aosoa, n_part);
  auto const ux = box_geo.length_inv()[0];
  auto const uy = box_geo.length_inv()[1];

  auto scatter_force = cs.get_scatter_force();
  kokkos_parallel_range_for<execution_space>(
      "ELC::far_field_forces", std::size_t{0u}, n_part, [&](std::size_t i) {
        auto const z = aosoa.position(i, 2);
        auto const q = aosoa.charge(i);
        Utils::Vector3d force{};
        for_each_entry(table, [&](FrequencyEntry const &entry, auto tag) {
          auto constexpr kind = decltype(tag)::value;
          auto const t = particle_terms<kind>(entry, scx, scy, i);
          auto const b = real_block<kind>(t, q, exp(entry.omega * z));
          auto const *g = coeffs.data() + entry.coeff_offset;
          if constexpr (kind == FreqKind::PQ) {
            // slots: 0 SSP, 1 SCP, 2 CSP, 3 CCP, 4 SSM, 5 SCM, 6 CSM, 7 CCM
            auto const pref_x =
                c_2pi * ux * static_cast<double>(entry.p) / entry.omega;
            auto const pref_y =
                c_2pi * uy * static_cast<double>(entry.q) / entry.omega;
            force[0] += pref_x *
                        (b[5] * g[3] + b[4] * g[2] - b[7] * g[1] - b[6] * g[0] +
                         b[1] * g[7] + b[0] * g[6] - b[3] * g[5] - b[2] * g[4]);
            force[1] += pref_y *
                        (b[6] * g[3] + b[4] * g[1] - b[7] * g[2] - b[5] * g[0] +
                         b[2] * g[7] + b[0] * g[5] - b[3] * g[6] - b[1] * g[4]);
            force[2] += (b[7] * g[3] + b[6] * g[2] + b[5] * g[1] + b[4] * g[0] -
                         b[3] * g[7] - b[2] * g[6] - b[1] * g[5] - b[0] * g[4]);
          } else {
            // slots: 0 SP, 1 CP, 2 SM, 3 CM
            auto constexpr axis = (kind == FreqKind::P) ? 0 : 1;
            force[axis] +=
                b[2] * g[1] - b[3] * g[0] + b[0] * g[3] - b[1] * g[2];
            force[2] += b[3] * g[1] + b[2] * g[0] - b[1] * g[3] - b[0] * g[2];
          }
        });
        auto access = scatter_force.access();
        access(i, 0) += force[0];
        access(i, 1) += force[1];
        access(i, 2) += force[2];
      });
}

double ElectrostaticLayerCorrection::calc_energy() const {
  auto const &system = get_system();
  auto const &box_geo = *system.box_geo;
  auto const &cs = *system.cell_structure;
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  auto energy = dipole_energy() + z_energy();

  auto const table = make_frequency_table(elc, prefactor, box_geo);
  auto const [n_freq_x, n_freq_y] = n_frequencies(box_geo, elc.far_cut);
  auto const scx =
      make_sc_cache(aosoa, n_part, 0, n_freq_x, box_geo.length_inv()[0]);
  auto const scy =
      make_sc_cache(aosoa, n_part, 1, n_freq_y, box_geo.length_inv()[1]);
  auto const coeffs =
      far_field_coefficients(elc, table, scx, scy, aosoa, n_part);

  energy += ELC::parallel_array_sum<1>(
      "ELC::far_field_energy", n_part,
      [&](std::size_t i, std::array<double, 1> &acc) {
        auto const z = aosoa.position(i, 2);
        auto const q = aosoa.charge(i);
        for_each_entry(table, [&](FrequencyEntry const &entry, auto tag) {
          auto constexpr kind = decltype(tag)::value;
          auto constexpr half = n_trig<kind>;
          auto const t = particle_terms<kind>(entry, scx, scy, i);
          auto const b = real_block<kind>(t, q, exp(entry.omega * z));
          auto const *g = coeffs.data() + entry.coeff_offset;
          // pair every P slot with its M slot and vice versa
          auto sum = 0.;
          for (std::size_t m = 0; m < half; ++m) {
            sum += b[half + m] * g[m] + b[m] * g[half + m];
          }
          acc[0] += sum / entry.omega;
        });
      })[0];
  /* we count both i<->j and j<->i, so return just half of it */
  return 0.5 * energy;
}

double ElectrostaticLayerCorrection::tune_far_cut() const {
  // Largest reasonable cutoff for far formula
  auto constexpr maximal_far_cut = 50.;
  auto const &box_geo = *get_system().box_geo;
  auto const box_l_x_inv = box_geo.length_inv()[0];
  auto const box_l_y_inv = box_geo.length_inv()[1];
  auto const min_inv_boxl = std::min(box_l_x_inv, box_l_y_inv);
  auto const box_l_z = box_geo.length()[2];
  auto const h = elc.box_h;
  // adjust lz according to dielectric layer method
  auto const lz = (elc.dielectric_contrast_on) ? h + elc.space_layer : box_l_z;

  auto tuned_far_cut = min_inv_boxl;
  double err;
  do {
    // following equation 18 in arnold02d
    auto const pref = 2. * std::numbers::pi * tuned_far_cut;
    auto const sum = pref + 2. * (box_l_x_inv + box_l_y_inv);
    auto const den = expm1(pref * lz);
    auto const num1 = exp(pref * h);
    auto const num2 = 1. / num1; // exp(-pref * h);

    err = 0.5 / den *
          (num1 / (lz - h) * (sum + 1. / (lz - h)) +
           num2 / (lz + h) * (sum + 1. / (lz + h)));

    tuned_far_cut += min_inv_boxl;
  } while (err > elc.maxPWerror and tuned_far_cut < maximal_far_cut);
  if (tuned_far_cut >= maximal_far_cut) {
    throw std::runtime_error("ELC tuning failed: maxPWerror too small");
  }
  return tuned_far_cut - min_inv_boxl;
}

void ElectrostaticLayerCorrection::sanity_checks_periodicity() const {
  auto const &box_geo = *get_system().box_geo;
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    throw std::runtime_error("ELC: requires periodicity (True, True, True)");
  }
}

void ElectrostaticLayerCorrection::adapt_solver() {
  std::visit(
      [this](auto &solver) {
        set_prefactor(solver->prefactor);
        solver->adapt_epsilon_elc();
        assert(solver->p3m_params.epsilon == P3M_EPSILON_METALLIC);
      },
      base_solver);
}

void ElectrostaticLayerCorrection::recalc_box_h() {
  m_box_geo = get_system().box_geo.get();
  auto const box_z = m_box_geo->length()[2];
  auto const new_box_h = box_z - elc.gap_size;
  if (new_box_h < 0.) {
    throw std::runtime_error("ELC gap size (" + std::to_string(elc.gap_size) +
                             ") larger than box length in z-direction (" +
                             std::to_string(box_z) + ")");
  }
  elc.box_h = new_box_h;
}

void ElectrostaticLayerCorrection::recalc_space_layer() {
  if (elc.dielectric_contrast_on) {
    auto const p3m_r_cut = std::visit(
        [](auto &solver) { return solver->p3m_params.r_cut; }, base_solver);
    // recalculate the space layer size:
    // 1. set the space_layer to be 1/3 of the gap size, so that box = layer
    elc.space_layer = (1. / 3.) * elc.gap_size;
    // 2. but make sure we don't overlap with the near-field formula
    auto const free_space = elc.gap_size - p3m_r_cut;
    // 3. and make sure the space layer is not bigger than half the actual
    // simulation box, to avoid overlaps
    auto const half_box_h = elc.box_h / 2.;
    auto const max_space_layer = std::min(free_space, half_box_h);
    if (elc.space_layer > max_space_layer) {
      if (max_space_layer <= 0.) {
        throw std::runtime_error("P3M real-space cutoff too large for ELC w/ "
                                 "dielectric contrast");
      }
      elc.space_layer = max_space_layer;
    }
    elc.space_box = elc.gap_size - 2. * elc.space_layer;
  }
}

elc_data::elc_data(double maxPWerror, double gap_size, double far_cut,
                   bool neutralize, double delta_top, double delta_bot,
                   bool with_const_pot, double potential_diff)
    : maxPWerror{maxPWerror}, gap_size{gap_size}, box_h{-1.}, far_cut{far_cut},
      far_cut2{-1.}, far_calculated{far_cut == -1.},
      dielectric_contrast_on{delta_top != 0. or delta_bot != 0.},
      const_pot{with_const_pot and dielectric_contrast_on},
      neutralize{neutralize and !dielectric_contrast_on},
      delta_mid_top{std::clamp(delta_top, -1., +1.)},
      delta_mid_bot{std::clamp(delta_bot, -1., +1.)},
      pot_diff{(with_const_pot) ? potential_diff : 0.},
      // initial setup of parameters, may change later when P3M is finally tuned
      // set the space_layer to be 1/3 of the gap size, so that box = layer
      space_layer{(dielectric_contrast_on) ? gap_size / 3. : 0.},
      space_box{gap_size - ((dielectric_contrast_on) ? 2. * space_layer : 0.)} {

  auto const delta_range = 1. + std::sqrt(round_error_prec);
  if (far_cut <= 0. and not far_calculated) {
    throw std::domain_error("Parameter 'far_cut' must be > 0");
  }
  if (maxPWerror <= 0.) {
    throw std::domain_error("Parameter 'maxPWerror' must be > 0");
  }
  if (gap_size <= 0.) {
    throw std::domain_error("Parameter 'gap_size' must be > 0");
  }
  if (potential_diff != 0. and not with_const_pot) {
    throw std::invalid_argument(
        "Parameter 'const_pot' must be True when 'pot_diff' is non-zero");
  }
  if (with_const_pot and not dielectric_contrast_on) {
    throw std::invalid_argument(
        "Parameter 'const_pot' requires a dielectric contrast; set "
        "'delta_mid_top' and 'delta_mid_bot' (use -1 for metallic walls)");
  }
  if (delta_top < -delta_range or delta_top > delta_range) {
    throw std::domain_error(
        "Parameter 'delta_mid_top' must be >= -1 and <= +1");
  }
  if (delta_bot < -delta_range or delta_bot > delta_range) {
    throw std::domain_error(
        "Parameter 'delta_mid_bot' must be >= -1 and <= +1");
  }
  /* Dielectric contrasts: the deltas should be either both -1 or both +1 when
   * no constant potential difference is applied. The case of two non-metallic
   * parallel boundaries can only be treated with a constant potential. */
  if (dielectric_contrast_on and not const_pot and
      (std::fabs(1. - delta_mid_top * delta_mid_bot) < round_error_prec)) {
    throw std::domain_error("ELC with two parallel metallic boundaries "
                            "requires the const_pot option");
  }
}

ElectrostaticLayerCorrection::ElectrostaticLayerCorrection(
    elc_data &&parameters, BaseSolver &&solver)
    : elc{parameters}, base_solver{solver} {
  // The P3M-GPU ignores images charges, disabled for now
  if (elc.dielectric_contrast_on) {
    auto const on_gpu =
        std::visit([](auto const &solver_ptr) { return solver_ptr->is_gpu(); },
                   base_solver);
    if (on_gpu) {
      throw std::runtime_error(
          "ELC with a dielectric contrast is not supported by the GPU "
          "variant of P3M");
    }
  }
  adapt_solver();
}

template <ChargeProtocol protocol>
void charge_assign(elc_data const &elc, CoulombP3M &solver) {
  auto constexpr assign_real =
      protocol == ChargeProtocol::BOTH or protocol == ChargeProtocol::REAL;
  auto constexpr assign_images =
      protocol == ChargeProtocol::BOTH or protocol == ChargeProtocol::IMAGE;
  solver.charge_assign_elc(
      {elc.space_layer, elc.box_h, elc.delta_mid_bot, elc.delta_mid_top},
      assign_real, assign_images);
}

template <ChargeProtocol protocol>
void modify_p3m_sums(elc_data const &elc, CoulombP3M &solver, auto const &cs) {
  auto const &aosoa = cs.get_aosoa();
  auto const n_part = cs.count_local_particles();
  // [0] number of charges, [1] sum q^2, [2] sum q
  auto const accumulate = [&elc, &aosoa](std::size_t i,
                                         std::array<double, 3> &acc) {
    auto const p_q = aosoa.charge(i);
    if (p_q == 0.) {
      return;
    }
    auto const p_z = aosoa.position(i, 2);
    if (protocol == ChargeProtocol::BOTH or protocol == ChargeProtocol::REAL) {
      acc[0] += 1.;
      acc[1] += Utils::sqr(p_q);
      acc[2] += p_q;
    }
    if (protocol == ChargeProtocol::BOTH or protocol == ChargeProtocol::IMAGE) {
      if (p_z < elc.space_layer) {
        acc[0] += 1.;
        acc[1] += Utils::sqr(elc.delta_mid_bot * p_q);
        acc[2] += elc.delta_mid_bot * p_q;
      }
      if (p_z > (elc.box_h - elc.space_layer)) {
        acc[0] += 1.;
        acc[1] += Utils::sqr(elc.delta_mid_top * p_q);
        acc[2] += elc.delta_mid_top * p_q;
      }
    }
  };
  // These sums feed P3M energies whose BOTH-minus-IMAGE difference cancels
  // strongly, so thread-order rounding would be amplified. The loop is O(N)
  // and cheap, therefore it stays serial.
  std::array<double, 3> local{};
  for (std::size_t i = 0; i < n_part; ++i) {
    accumulate(i, local);
  }
  std::array<double, 3> global{};
  boost::mpi::all_reduce(comm_cart, local.data(), 3, global.data(),
                         std::plus<>());
  auto const global_n = static_cast<std::size_t>(std::llround(global[0]));
  solver.count_charged_particles_elc(global_n, global[1],
                                     Utils::sqr(global[2]));
}

double ElectrostaticLayerCorrection::long_range_energy() const {
  auto const &system = get_system();
  auto const energy = std::visit(
      [this, &system](auto const &solver_ptr) {
        auto &solver = *solver_ptr;
        auto const &cs = *system.cell_structure;
        auto const &box_geo = *system.box_geo;

        // assign the original charges (they may not have been assigned yet)
        solver.charge_assign();

        if (!elc.dielectric_contrast_on) {
          return solver.long_range_energy();
        }

        auto energy = 0.;
        energy += 0.5 * solver.long_range_energy();
        energy += 0.5 * elc.dielectric_layers_self_energy(solver, box_geo, cs);

        // assign both original and image charges
        charge_assign<ChargeProtocol::BOTH>(elc, solver);
        modify_p3m_sums<ChargeProtocol::BOTH>(elc, solver, cs);
        energy += 0.5 * solver.long_range_energy();

        // assign only the image charges now
        charge_assign<ChargeProtocol::IMAGE>(elc, solver);
        modify_p3m_sums<ChargeProtocol::IMAGE>(elc, solver, cs);
        energy -= 0.5 * solver.long_range_energy();

        // restore modified sums
        modify_p3m_sums<ChargeProtocol::REAL>(elc, solver, cs);

        return energy;
      },
      base_solver);
  return energy + calc_energy();
}

void ElectrostaticLayerCorrection::add_long_range_forces() const {
  auto const &system = get_system();
  std::visit(
      [this, &system](auto const &solver_ptr) {
        auto &cs = *system.cell_structure;
        auto &solver = *solver_ptr;
        if (elc.dielectric_contrast_on) {
          auto const &box_geo = *system.box_geo;
          modify_p3m_sums<ChargeProtocol::BOTH>(elc, solver, cs);
          charge_assign<ChargeProtocol::BOTH>(elc, solver);
          elc.dielectric_layers_self_forces(solver, box_geo, cs);
        } else {
          solver.charge_assign();
        }
        solver.add_long_range_forces();
        if (elc.dielectric_contrast_on) {
          modify_p3m_sums<ChargeProtocol::REAL>(elc, solver, cs);
        }
      },
      base_solver);
  add_force();
}

#endif // ESPRESSO_P3M
