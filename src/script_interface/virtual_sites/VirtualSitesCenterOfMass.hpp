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

#ifndef SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_CENTER_OF_MASS_HPP
#define SCRIPT_INTERFACE_VIRTUAL_SITES_VIRTUAL_SITES_CENTER_OF_MASS_HPP

#include "config/config.hpp"

#ifdef VIRTUAL_SITES_CENTER_OF_MASS

#include "VirtualSites.hpp"

#include "core/virtual_sites/VirtualSitesCenterOfMass.hpp"

#include <memory>

namespace ScriptInterface {
namespace VirtualSites {

class VirtualSitesCenterOfMass : public VirtualSites {
public:
  VirtualSitesCenterOfMass() : m_virtual_sites(std::make_shared<::VirtualSitesCenterOfMass>()) {
    // add_parameters({
    //   {"mid_for_vs", AutoParameter::read_only,
    //       [this]() { return Variant( m_virtual_sites->get_mid_for_vs()); }}
    // });
  }

  /** Vs implementation we are wrapping */
  std::shared_ptr<::VirtualSites> virtual_sites() override {
    return m_virtual_sites;
  } // TO DO: MAKE THAT THIS CONSTRUCTOR TAKES AN UNDORDERED_MAP TO
    // INSTANTIATE THE VIRTUAL_SITES_CENTER_OF_MASS

private:
  std::shared_ptr<::VirtualSitesCenterOfMass> m_virtual_sites;
};

} /* namespace VirtualSites */
} /* namespace ScriptInterface */
#endif // VIRTUAL_SITES_CENTER_OF_MASS
#endif
