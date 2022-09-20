/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_COULOMB_HPP

#include "config/config.hpp"

#ifdef ELECTROSTATICS

#include "actor/traits.hpp"

#include "electrostatics/debye_hueckel.hpp"
#include "electrostatics/elc.hpp"
#include "electrostatics/icc.hpp"
#include "electrostatics/mmm1d.hpp"
#include "electrostatics/mmm1d_gpu.hpp"
#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m_gpu.hpp"
#include "electrostatics/reaction_field.hpp"
#include "electrostatics/scafacos.hpp"

#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <cstddef>
#include <functional>
#include <memory>
#include <type_traits>

using ElectrostaticsActor =
    boost::variant<std::shared_ptr<DebyeHueckel>,
#ifdef P3M
                   std::shared_ptr<CoulombP3M>,
#ifdef CUDA
                   std::shared_ptr<CoulombP3MGPU>,
#endif // CUDA
                   std::shared_ptr<ElectrostaticLayerCorrection>,
#endif // P3M
                   std::shared_ptr<CoulombMMM1D>,
#ifdef MMM1D_GPU
                   std::shared_ptr<CoulombMMM1DGpu>,
#endif // MMM1D_GPU
#ifdef SCAFACOS
                   std::shared_ptr<CoulombScafacos>,
#endif // SCAFACOS
                   std::shared_ptr<ReactionField>>;

using ElectrostaticsExtension = boost::variant<std::shared_ptr<ICCStar>>;

extern boost::optional<ElectrostaticsActor> electrostatics_actor;
extern boost::optional<ElectrostaticsExtension> electrostatics_extension;

/** Get the electrostatics prefactor. */
struct GetCoulombPrefactor : public boost::static_visitor<double> {
  template <typename T>
  double operator()(std::shared_ptr<T> const &actor) const {
    return actor->prefactor;
  }
};

namespace Coulomb {

namespace traits {

#ifdef P3M
/** @brief Whether an actor can be adapted by ELC. */
template <typename T>
using elc_adaptable =
    std::is_convertible<std::shared_ptr<T>,
                        ElectrostaticLayerCorrection::BaseSolver>;
#endif // P3M

/** @brief Whether an actor is a solver. */
template <typename T>
using is_solver = std::is_convertible<std::shared_ptr<T>, ElectrostaticsActor>;
/** @brief Whether an actor is an extension. */
template <typename T>
using is_extension =
    std::is_convertible<std::shared_ptr<T>, ElectrostaticsExtension>;

/** @brief The electrostatic method supports pressure calculation. */
template <class T> struct has_pressure : std::true_type {};
#ifdef P3M
template <>
struct has_pressure<ElectrostaticLayerCorrection> : std::false_type {};
#endif // P3M
#ifdef MMM1D_GPU
template <> struct has_pressure<CoulombMMM1DGpu> : std::false_type {};
#endif // MMM1D_GPU
#ifdef SCAFACOS
template <> struct has_pressure<CoulombScafacos> : std::false_type {};
#endif // SCAFACOS
template <> struct has_pressure<CoulombMMM1D> : std::false_type {};

} // namespace traits

/** @brief Check if the system is charge-neutral. */
void check_charge_neutrality(double relative_tolerance);

Utils::Vector9d calc_pressure_long_range(ParticleRange const &particles);

void sanity_checks();
double cutoff();

void on_observable_calc();
void on_coulomb_change();
void on_boxl_change();
void on_node_grid_change();
void on_periodicity_change();
void on_cell_structure_change();

void calc_long_range_force(ParticleRange const &particles);
double calc_energy_long_range(ParticleRange const &particles);

namespace detail {
bool flag_all_reduce(bool flag);
} // namespace detail

} // namespace Coulomb
#else // ELECTROSTATICS
#include <cstddef>
namespace Coulomb {
constexpr std::size_t pressure_n() { return 0; }
constexpr std::size_t energy_n() { return 0; }
} // namespace Coulomb
#endif // ELECTROSTATICS
#endif
