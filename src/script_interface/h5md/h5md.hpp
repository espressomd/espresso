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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_H5MD_H5MD_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_H5MD_H5MD_HPP

#include "config/config.hpp"

#ifdef H5MD

#include "io/writer/h5md_core.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Writer {
class H5md : public AutoParameters<H5md> {
public:
  H5md() {
    add_parameters(
        {{"file_path", m_h5md, &::Writer::H5md::File::file_path},
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

private:
  Variant do_call_method(const std::string &name,
                         const VariantMap &parameters) override;

  void do_construct(VariantMap const &params) override {
    m_output_fields = get_value<std::vector<std::string>>(params, "fields");
    m_h5md = make_shared_from_args<::Writer::H5md::File, std::string,
                                   std::string, std::vector<std::string>,
                                   std::string, std::string, std::string,
                                   std::string, std::string, std::string>(
        params, "file_path", "script_path", "fields", "mass_unit",
        "length_unit", "time_unit", "force_unit", "velocity_unit",
        "charge_unit");
  }

  std::shared_ptr<::Writer::H5md::File> m_h5md;
  std::vector<std::string> m_output_fields;
};

} // namespace Writer
} // namespace ScriptInterface

#endif // H5MD
#endif
