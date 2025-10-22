/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include "AccumulatorBase.hpp"
#include "observables/Observable.hpp"
#include "system/System.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Accumulators {

/**
 * @brief Record the time distances are below a threshold.
 *
 * The contact time is defined as @f$ \tau = t_f - t_0 @f$ where
 * @f$ t_f @f$ is the last measured time at which a pairwise distance was below
 * the threshold and @f$ t_0 @f$ is the first measured time at which the same
 * distance was below the threshold.
 */
class ContactTimes : public AccumulatorBase {
public:
  ContactTimes(::System::System const *system, int delta_N,
               double contact_threshold,
               std::shared_ptr<Observables::Observable> obs)
      : AccumulatorBase(system, delta_N), m_obs(std::move(obs)),
        m_contact_threshold(contact_threshold) {
    if (contact_threshold < 0.) {
      throw std::domain_error("Attribute 'contact_threshold' must be >= 0");
    }
  }

  void update(boost::mpi::communicator const &comm) override;
  std::string get_internal_state() const final;
  void set_internal_state(std::string const &) final;

  auto const &contact_times() const { return m_data; }
  std::vector<std::size_t> shape() const override { return {m_data.size()}; }
  void clear() { m_data.clear(); }
  double get_contact_threshold() const { return m_contact_threshold; }

private:
  std::shared_ptr<Observables::Observable> m_obs;
  std::vector<double> m_data;
  double m_contact_threshold;
  std::vector<bool> contacts;
  std::vector<double> first_contact_times;
};

} // namespace Accumulators
