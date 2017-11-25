/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project

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

#include "Observable.hpp"
#include "integrate.hpp"
#include "partCfg_global.hpp"
#include <fstream>
#include <string>

namespace Observables {

Observable::Observable() : m_ofile(), m_filename(), m_binary(false) {}

void Observable::set_filename(std::string const &filename, bool binary) {
  if (!filename.empty()) {
    m_filename = filename;
    m_binary = binary;

    auto mode = std::ios_base::app;
    if (m_binary)
      mode |= std::ios_base::binary;

    if (m_ofile.is_open()) {
      m_ofile.close();
      m_ofile.clear(); // clear flags
    }

    m_ofile.open(m_filename, mode);
  }
}

bool Observable::writable() const { return m_ofile.is_open(); }

void Observable::write() {
  if (writable()) {
    do_write();
  }
}

void Observable::do_write() {
  m_ofile << sim_time;
  for (auto p : operator()(partCfg()))
    m_ofile << " " << p;
  m_ofile << std::endl;
}

} // Namespace Observables
