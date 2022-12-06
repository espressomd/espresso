/*
 * Copyright (C) 2022 The ESPResSo project
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

#pragma once

#include "walberla_bridge/electrokinetics/PoissonSolver/PoissonSolver.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <vector>

template <class EKSpecies> class EKContainer {
  using container_type = std::vector<std::shared_ptr<EKSpecies>>;

public:
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

private:
  container_type m_ekcontainer;
  double m_tau{};

  // TODO: this could be moved to the Scriptinterface
  std::shared_ptr<walberla::PoissonSolver> m_poissonsolver;

public:
  void add(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) ==
           m_ekcontainer.end());

    m_ekcontainer.emplace_back(c);
  }
  void remove(std::shared_ptr<EKSpecies> const &c) {
    assert(std::find(m_ekcontainer.begin(), m_ekcontainer.end(), c) !=
           m_ekcontainer.end());
    m_ekcontainer.erase(
        std::remove(m_ekcontainer.begin(), m_ekcontainer.end(), c),
        m_ekcontainer.end());
  }

  iterator begin() noexcept { return m_ekcontainer.begin(); }
  iterator end() noexcept { return m_ekcontainer.end(); }
  const_iterator begin() const noexcept { return m_ekcontainer.begin(); }
  const_iterator end() const noexcept { return m_ekcontainer.end(); }
  [[nodiscard]] bool empty() const noexcept { return m_ekcontainer.empty(); }

  void set_poissonsolver(
      std::shared_ptr<walberla::PoissonSolver> const &solver) noexcept {
    m_poissonsolver = solver;
  }

  [[nodiscard]] bool is_poissonsolver_set() const noexcept {
    return m_poissonsolver != nullptr;
  }

  [[nodiscard]] double get_tau() const noexcept { return m_tau; }
  void set_tau(double tau) noexcept { m_tau = tau; }

  void reset_charge() const { m_poissonsolver->reset_charge_field(); }
  void add_charge(std::size_t const id, double valency,
                  bool is_double_precision) const {
    m_poissonsolver->add_charge_to_field(id, valency, is_double_precision);
  }

  void solve_poisson() const { m_poissonsolver->solve(); }

  [[nodiscard]] std::size_t get_potential_field_id() const {
    return m_poissonsolver->get_potential_field_id();
  }
};
