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

#include "config/config.hpp"

#ifdef ESPRESSO_H5MD

#include "h5md.hpp"

#include "core/cell_system/CellStructure.hpp"
#include "core/communication.hpp"
#include "core/io/writer/h5md_core.hpp"
#include "core/system/System.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Writer {

H5md::H5md() {
  add_parameters(
      {{"file_path", m_h5md, &::Writer::H5md::File::file_path},
       {"chunk_size", m_h5md, &::Writer::H5md::File::chunk_size},
       {"script_path", m_h5md, &::Writer::H5md::File::script_path},
       {"fields", AutoParameter::read_only,
        [this]() { return make_vector_of_variants(m_output_fields); }},
       {"mass_unit", m_h5md, &::Writer::H5md::File::mass_unit},
       {"length_unit", m_h5md, &::Writer::H5md::File::length_unit},
       {"time_unit", m_h5md, &::Writer::H5md::File::time_unit},
       {"force_unit", m_h5md, &::Writer::H5md::File::force_unit},
       {"velocity_unit", m_h5md, &::Writer::H5md::File::velocity_unit},
       {"charge_unit", m_h5md, &::Writer::H5md::File::charge_unit}});
};

void H5md::do_construct(VariantMap const &params) {
  m_output_fields = get_value<std::vector<std::string>>(params, "fields");
  m_h5md =
      make_shared_from_args<::Writer::H5md::File, std::filesystem::path,
                            std::filesystem::path, std::vector<std::string>,
                            std::string, std::string, std::string, std::string,
                            std::string, std::string, int>(
          params, "file_path", "script_path", "fields", "mass_unit",
          "length_unit", "time_unit", "force_unit", "velocity_unit",
          "charge_unit", "chunk_size");
  // MPI communicator is needed to close parallel file handles
  m_mpi_env_lock = ::communication_environment->get_mpi_env();
  is_file_closed = false;
}

H5md::~H5md() {
  m_h5md.reset();
  assert(m_h5md.use_count() == 0u);
  m_mpi_env_lock.reset();
}

void H5md::sanity_check_is_file_closed(std::string const &name) const {
  if (is_file_closed) {
    if (context()->is_head_node()) {
      throw std::runtime_error("cannot call '" + name + "' on a closed file");
    }
    throw Exception("");
  }
}

Variant H5md::do_call_method(std::string const &name, VariantMap const &) {
  if (name == "write") {
    sanity_check_is_file_closed(name);
    auto const &system = ::System::get_system();
    // Forces/torques live in the ParticleStore columns, which are only rebuilt
    // lazily. A write triggered before the first integrator step (or after a
    // particle was added) would otherwise read stale/detached rows and crash;
    // sync the store first.
    system.cell_structure->ensure_particle_store_synchronized();
    auto const particles = system.cell_structure->local_particles();
    auto const sim_time = system.get_sim_time();
    auto const time_step = system.get_time_step();
    auto const n_steps = static_cast<int>(std::round(sim_time / time_step));
    m_h5md->write(particles, sim_time, n_steps, *system.box_geo);
  } else if (name == "flush") {
    sanity_check_is_file_closed(name);
    m_h5md->flush();
  } else if (name == "close") {
    sanity_check_is_file_closed(name);
    m_h5md->close();
    is_file_closed = true;
  } else if (name == "valid_fields") {
    return make_vector_of_variants(m_h5md->valid_fields());
  }
  return {};
}

} // namespace Writer
} // namespace ScriptInterface

#endif // ESPRESSO_H5MD
