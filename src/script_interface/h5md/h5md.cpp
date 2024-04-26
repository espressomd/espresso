/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifdef H5MD

#include "h5md.hpp"

#include "cell_system/CellStructure.hpp"
#include "core/system/System.hpp"

#include <cmath>
#include <string>

namespace ScriptInterface {
namespace Writer {

Variant H5md::do_call_method(const std::string &name,
                             const VariantMap &parameters) {
  if (name == "write") {
    auto const &system = ::System::get_system();
    auto const particles = system.cell_structure->local_particles();
    auto const sim_time = system.get_sim_time();
    auto const time_step = system.get_time_step();
    auto const n_steps = static_cast<int>(std::round(sim_time / time_step));
    m_h5md->write(particles, sim_time, n_steps, *system.box_geo);
  } else if (name == "flush") {
    m_h5md->flush();
  } else if (name == "close") {
    m_h5md->close();
  } else if (name == "valid_fields") {
    return make_vector_of_variants(m_h5md->valid_fields());
  }
  return {};
}

} // namespace Writer
} // namespace ScriptInterface

#endif // H5MD
