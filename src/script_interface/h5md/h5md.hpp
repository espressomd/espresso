/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#pragma once

#include "config/config.hpp"

#ifdef ESPRESSO_H5MD

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <boost/mpi/environment.hpp>

#include <memory>
#include <string>
#include <vector>

namespace Writer::H5md {
class File;
} // namespace Writer::H5md

namespace ScriptInterface {
namespace Writer {

class H5md : public AutoParameters<H5md> {
public:
  H5md();

  Variant do_call_method(const std::string &name,
                         const VariantMap &parameters) override;

  void do_construct(VariantMap const &params) override;

  ~H5md() override;

private:
  std::shared_ptr<boost::mpi::environment> m_mpi_env_lock;
  std::shared_ptr<::Writer::H5md::File> m_h5md;
  std::vector<std::string> m_output_fields;
};

} // namespace Writer
} // namespace ScriptInterface

#endif // ESPRESSO_H5MD
