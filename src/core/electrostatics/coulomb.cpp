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

#include "config/config.hpp"

#include "electrostatics/solver.hpp"

#ifdef ELECTROSTATICS

#include "electrostatics/coulomb.hpp"

#include "ParticleRange.hpp"
#include "actor/visit_try_catch.hpp"
#include "actor/visitors.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "electrostatics/icc.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/checks/charge_neutrality.hpp>
#include <utils/constants.hpp>
#include <utils/demangle.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>
#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/gather.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

namespace Coulomb {

Solver::Solver() {
  impl = std::make_unique<Implementation>();
  reinit_on_observable_calc = false;
}

Solver const &get_coulomb() { return System::get_system().coulomb; }

void Solver::sanity_checks() const {
  if (impl->solver) {
    std::visit([](auto const &ptr) { ptr->sanity_checks(); }, *impl->solver);
  }
}

void Solver::on_coulomb_change() {
  reinit_on_observable_calc = true;
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->init(); }, *impl->solver);
  }
}

void Solver::on_boxl_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_boxl_change(); }, *impl->solver);
  }
}

void Solver::on_node_grid_change() {
  if (impl->solver) {
    std::visit([](auto &ptr) { ptr->on_node_grid_change(); }, *impl->solver);
  }
}

void Solver::on_periodicity_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_periodicity_change(); },
                    *impl->solver);
  }
}

void Solver::on_cell_structure_change() {
  if (impl->solver) {
    visit_try_catch([](auto &ptr) { ptr->on_cell_structure_change(); },
                    *impl->solver);
  }
}

struct LongRangePressure {
  explicit LongRangePressure(ParticleRange const &particles)
      : m_particles{particles} {}

#ifdef P3M
  auto operator()(std::shared_ptr<CoulombP3M> const &actor) const {
    return actor->long_range_pressure(m_particles);
  }
#endif // P3M

  auto operator()(std::shared_ptr<DebyeHueckel> const &actor) const {
    return Utils::Vector9d{};
  }

  auto operator()(std::shared_ptr<ReactionField> const &actor) const {
    return Utils::Vector9d{};
  }

  template <typename T,
            std::enable_if_t<!traits::has_pressure<T>::value> * = nullptr>
  auto operator()(std::shared_ptr<T> const &) const {
    runtimeWarningMsg() << "Pressure calculation not implemented by "
                        << "electrostatics method " << Utils::demangle<T>();
    return Utils::Vector9d{};
  }

private:
  ParticleRange const &m_particles;
};

Utils::Vector9d
Solver::calc_pressure_long_range(ParticleRange const &particles) const {
  if (impl->solver) {
    return std::visit(LongRangePressure(particles), *impl->solver);
  }
  return {};
}

struct ShortRangeCutoff {
#ifdef P3M
  auto operator()(std::shared_ptr<CoulombP3M> const &actor) const {
    return actor->p3m.params.r_cut;
  }
  auto
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    return std::max(actor->elc.space_layer,
                    std::visit(*this, actor->base_solver));
  }
#endif // P3M
#ifdef MMM1D_GPU
  auto operator()(std::shared_ptr<CoulombMMM1DGpu> const &actor) const {
    return std::numeric_limits<double>::infinity();
  }
#endif // MMM1D_GPU
  auto operator()(std::shared_ptr<CoulombMMM1D> const &actor) const {
    return std::numeric_limits<double>::infinity();
  }
#ifdef SCAFACOS
  auto operator()(std::shared_ptr<CoulombScafacos> const &actor) const {
    return actor->get_r_cut();
  }
#endif // SCAFACOS
  auto operator()(std::shared_ptr<ReactionField> const &actor) const {
    return actor->r_cut;
  }
  auto operator()(std::shared_ptr<DebyeHueckel> const &actor) const {
    return actor->r_cut;
  }
};

double Solver::cutoff() const {
  if (impl->solver) {
    return std::visit(ShortRangeCutoff(), *impl->solver);
  }
  return -1.0;
}

struct EventOnObservableCalc {
  template <typename T> void operator()(std::shared_ptr<T> const &) const {}

#ifdef P3M
  void operator()(std::shared_ptr<CoulombP3M> const &actor) const {
    actor->count_charged_particles();
  }
  void
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    std::visit(*this, actor->base_solver);
  }
#endif // P3M
};

void Solver::on_observable_calc() {
  if (reinit_on_observable_calc) {
    if (impl->solver) {
      std::visit(EventOnObservableCalc(), *impl->solver);
    }
    reinit_on_observable_calc = false;
  }
}

struct LongRangeForce {
  explicit LongRangeForce(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef P3M
  void operator()(std::shared_ptr<CoulombP3M> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#ifdef CUDA
  void operator()(std::shared_ptr<CoulombP3MGPU> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#endif // CUDA
  void
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    actor->add_long_range_forces(m_particles);
  }
#endif // P3M
#ifdef MMM1D_GPU
  void operator()(std::shared_ptr<CoulombMMM1DGpu> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
#ifdef SCAFACOS
  void operator()(std::shared_ptr<CoulombScafacos> const &actor) const {
    actor->add_long_range_forces();
  }
#endif
  /* Several algorithms only provide near-field kernels */
  void operator()(std::shared_ptr<CoulombMMM1D> const &) const {}
  void operator()(std::shared_ptr<DebyeHueckel> const &) const {}
  void operator()(std::shared_ptr<ReactionField> const &) const {}

private:
  ParticleRange const &m_particles;
};

struct LongRangeEnergy {
  explicit LongRangeEnergy(ParticleRange const &particles)
      : m_particles(particles) {}

#ifdef P3M
  auto operator()(std::shared_ptr<CoulombP3M> const &actor) const {
    return actor->long_range_energy(m_particles);
  }
  auto
  operator()(std::shared_ptr<ElectrostaticLayerCorrection> const &actor) const {
    return actor->long_range_energy(m_particles);
  }
#endif // P3M
#ifdef MMM1D_GPU
  auto operator()(std::shared_ptr<CoulombMMM1DGpu> const &actor) const {
    actor->add_long_range_energy();
    return 0.;
  }
#endif // MMM1D_GPU
#ifdef SCAFACOS
  auto operator()(std::shared_ptr<CoulombScafacos> const &actor) const {
    return actor->long_range_energy();
  }
#endif
  /* Several algorithms only provide near-field kernels */
  auto operator()(std::shared_ptr<CoulombMMM1D> const &) const { return 0.; }
  auto operator()(std::shared_ptr<DebyeHueckel> const &) const { return 0.; }
  auto operator()(std::shared_ptr<ReactionField> const &) const { return 0.; }

private:
  ParticleRange const &m_particles;
};

void Solver::calc_long_range_force(ParticleRange const &particles) const {
  if (impl->solver) {
    std::visit(LongRangeForce(particles), *impl->solver);
  }
}

double Solver::calc_energy_long_range(ParticleRange const &particles) const {
  if (impl->solver) {
    return std::visit(LongRangeEnergy(particles), *impl->solver);
  }
  return 0.;
}

/** @brief Compute the net charge rescaled by the smallest non-zero charge. */
static auto calc_charge_excess_ratio(std::vector<double> const &charges) {
  using namespace boost::accumulators;
  using KahanSum = accumulator_set<double, features<tag::sum_kahan>>;

  KahanSum q_sum;
  auto q_min = std::numeric_limits<double>::infinity();

  for (auto const q : charges) {
    if (q != 0.) {
      q_sum(q);
      q_min = std::min(q_min, std::abs(q));
    }
  }

  return std::abs(sum_kahan(q_sum)) / q_min;
}

void check_charge_neutrality(System::System const &system,
                             double relative_tolerance) {
  // collect non-zero particle charges from all nodes
  std::vector<double> local_charges;
  for (auto const &p : system.cell_structure->local_particles()) {
    local_charges.push_back(p.q());
  }
  std::vector<std::vector<double>> node_charges;
  boost::mpi::gather(comm_cart, local_charges, node_charges, 0);

  // run Kahan sum on charges
  auto excess_ratio = 0.;
  if (this_node == 0) {
    auto charges = std::move(local_charges);
    for (auto it = ++node_charges.begin(); it != node_charges.end(); ++it) {
      charges.insert(charges.end(), it->begin(), it->end());
    }
    excess_ratio = calc_charge_excess_ratio(charges);
  }
  boost::mpi::broadcast(comm_cart, excess_ratio, 0);

  if (excess_ratio >= relative_tolerance) {
    std::ostringstream serializer;
    serializer << std::scientific << std::setprecision(4);
    serializer << excess_ratio;
    throw std::runtime_error(
        "The system is not charge neutral. Please neutralize the system "
        "before adding a new actor by adding the corresponding counterions "
        "to the system. Alternatively you can turn off the electroneutrality "
        "check by supplying check_neutrality=False when creating the actor. "
        "In this case you may be simulating a non-neutral system which will "
        "affect physical observables like e.g. the pressure, the chemical "
        "potentials of charged species or potential energies of the system. "
        "Since simulations of non charge neutral systems are special please "
        "make sure you know what you are doing. The relative charge excess "
        "was " +
        serializer.str());
  }
}

} // namespace Coulomb
#endif // ELECTROSTATICS
