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

#include "ContactTimes.hpp"

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/vector.hpp>

#include <sstream>
#include <string>

namespace Accumulators {
void ContactTimes::update(boost::mpi::communicator const &comm) {
  if (comm.rank() != 0) {
    return;
  }
  auto const &system = System::get_system();
  auto time = system.get_sim_time();
  auto dt = system.get_time_step();
  auto dists = m_obs->operator()(comm);

  // Initialize the bookkeeping vectors
  if (contacts.empty()) {
    contacts.resize(dists.size(), false);
  }

  if (first_contact_times.empty()) {
    first_contact_times.resize(dists.size(), time);
  }

  for (std::size_t i = 0u; i < dists.size(); ++i) {
    if (dists[i] < m_contact_threshold) { // distance below the threshold
      if (not contacts[i]) {              // but it was not before!
        contacts[i] = true;
        first_contact_times[i] = time;
      }
    } else {             // distance above the threshold
      if (contacts[i]) { // it was below before
        // calculate the total contact time
        auto contact_time = time - first_contact_times[i] - dt;
        if (contact_time < dt) {
          contact_time = 0.;
        }
        contacts[i] = false;
        m_data.emplace_back(contact_time);
      }
    }
  }
}

std::string ContactTimes::get_internal_state() const {
  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << m_data;

  return ss.str();
}

void ContactTimes::set_internal_state(std::string const &state) {
  namespace iostreams = boost::iostreams;
  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);
  boost::archive::binary_iarchive ia(ss);

  ia >> m_data;
  m_system = nullptr;
}
} // namespace Accumulators
