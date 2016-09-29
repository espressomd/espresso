/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef ESPRESSO_H5MD_HPP
#define ESPRESSO_H5MD_HPP

#include "ScriptInterface.hpp"

namespace ScriptInterface {
namespace Writer {
namespace H5md {

class H5md : public ScriptInterfaceBase {
    public:
        // Returns the name of the class
        const std::string name() const override { return
                "ScriptInterface::Writer::H5md::H5md"; }
        Variant call_method(std::string const &method,
                            VariantMap const &parameters);
};

} /* namespace H5md */
} /* namespace Writer */
} /* namespace Scriptinterface */


#endif //ESPRESSO_H5MD_HPP
