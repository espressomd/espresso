/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#ifndef __SCRIPT_INTERFACE_HPP
#define __SCRIPT_INTERFACE_HPP

#include <map>
#include <string>
#include <vector>

#include <boost/variant.hpp>

/**
 * @brief Possible types for parameters.
 */
typedef boost::variant<bool,
                       int,
                       double,
                       std::string,
                       std::vector<int>,
                       std::vector<double> > Variant;

#include "Parameter.hpp"

/**
 * @brief Tries to extract a value with the type of MEMBER_NAME from the Variant.
 *
 * This will fail at compile time if the type of MEMBER_NAME is not one of the
 * possible types of Variant, and at runtime if the current type of the variant
 * is not that of MEMBER_NAME.
 */
#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME) if(name == PARAMETER_NAME) {\
  MEMBER_NAME = boost::get<decltype(MEMBER_NAME)>(value);\
  }

namespace ScriptInterface {

/**
 * @brief Make a Variant from argument.
 * 
 * This is a convinience function, so that
 * rather involved constructors from
 * boost::variant are not needed in the
 * script interfaces.
 */
  template<typename T>
  Variant make_variant(const T& x) {
    return Variant(x);
  }

/**
 * @brief Base class for generic script interface.
 *
 * @TODO Add extensive documentation.
 *
 */
class ScriptInterfaceBase {
 public:
  /**
   * @brief Human-readable name of the object.
   *
   * @return Name of the object.
   */   
  virtual const std::string name() const = 0;

  /**
   * @brief get current parameters.
   * This is not a reference on purpose becaue
   * operator[] is not const in std::map.
   *
   * @return Parameters set in class.
   */
  virtual std::map<std::string, Variant> get_parameters() const = 0;

  /**
   * @brief Get requiered and optional parameters for class
   *
   * Get requiered and optional parameters for class.
   *
   * @return Expected parameters.
   */  
  virtual std::map<std::string, Parameter> all_parameters() const = 0;
  
  /**
   * @brief Get single parameter.
   *
   * @param name Name of the parameter
   * @return Value of parameter @param name.
   */
  virtual Variant get_parameter(const std::string &name) const {
    return get_parameters().at(name);
  }
  
  /**
   * @brief Set single parameter.
   *
   * @param name Name of the parameter
   * @param value Set parameter to this value.
   */  
  virtual void set_parameter(const std::string &name, const Variant &value) = 0;
  /**
   * @brief Set multiple parameters.
   *
   * The default implementation calls the implementation of set_parameter for every
   * element of the map.
   *
   * @param Paramters Parameters to set.
   */  
  virtual void set_parameters(const std::map<std::string, Variant> &parameters)
  {
    for (auto const& it : parameters) {
      set_parameter(it.first, it.second);
    }
  }
};

} /* namespace ScriptInterface */

#endif
