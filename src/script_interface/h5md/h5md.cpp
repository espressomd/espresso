/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "config.hpp"
#ifdef H5MD
#ifndef ESPRESSO_SCRIPTINTERFACE_H5MD_CPP
#define ESPRESSO_SCRIPTINTERFACE_H5MD_CPP
#include "cells.hpp"
#include "h5md.hpp"
#include "partCfg_global.hpp"

namespace ScriptInterface {
namespace Writer {
Variant H5mdScript::call_method(const std::string &name,
                                const VariantMap &parameters) {
  if (name == "init_file")
    m_h5md->InitFile();
  else if (name == "write")
    m_h5md->Write(m_h5md->what(), partCfg(), cell_structure.local_particles());
  else if (name == "flush")
    m_h5md->Flush();
  else if (name == "close")
    m_h5md->Close();
  return {};
}

} /* namespace Writer */
} // namespace ScriptInterface

#endif // ESPRESSO_H5MD_HPP
#endif // H5MD
