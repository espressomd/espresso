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

#ifndef SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP
#define SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP

#include <map>
#include <string>
#include <type_traits>
#include <vector>

#include <boost/variant.hpp>

#include "utils/NumeratedContainer.hpp"

#include "Parameter.hpp"
#include "core/Vector.hpp"

namespace Utils {
template <typename T> class ObjectId {
public:
  /* Assign an id on construction */
  ObjectId() : m_id(reg().add(std::weak_ptr<T>())) {}

  /* Remove id on destruction */
  virtual ~ObjectId() { reg().remove(m_id); }

  /**
   * @brief Get indentifier for this instance.
   */
  int id() const { return m_id; }
  /**
   * @brief get instance by id.
   */
  static std::weak_ptr<T> &get_instance(int id) { return reg()[id]; }

private:
  const int m_id;
  static Utils::NumeratedContainer<std::weak_ptr<T>> &reg() {
    static Utils::NumeratedContainer<std::weak_ptr<T>> m_reg;

    return m_reg;
  }
};
} /* namespace Utils */

namespace ScriptInterface {

struct OId {
  OId() : id(-1) {}
  explicit OId(int i) : id(i) {}

  int id;

  template <typename Archive>
  void serialize(Archive &ar, unsigned int /* version */) {
    ar &id;
  }
};

std::ostream &operator<<(std::ostream &out, OId const &oid);
/**
 * @brief Possible types for parameters.
 */
typedef boost::variant<bool, int, double, std::string, std::vector<int>,
                       std::vector<double>, Vector2d, Vector3d, OId>
    Variant;

/**
 * @brief Tries to extract a value with the type of MEMBER_NAME from the
 * Variant.
 *
 * This will fail at compile time if the type of MEMBER_NAME is not one of the
 * possible types of Variant, and at runtime if the current type of the variant
 * is not that of MEMBER_NAME.
 * remove_reference ensures that this also works with member access by reference
 * for example as returned by a function.
 */

#define SET_PARAMETER_HELPER(PARAMETER_NAME, MEMBER_NAME)                      \
  if (name == PARAMETER_NAME) {                                                \
    MEMBER_NAME =                                                              \
        boost::get<std::remove_reference<decltype(MEMBER_NAME)>::type>(value); \
  }

#define SET_PARAMETER_HELPER_VECTOR3D(PARAMETER_NAME, MEMBER_NAME)             \
  if (name == PARAMETER_NAME) {                                                \
    MEMBER_NAME = Vector3d(boost::get<std::vector<double>>(value).data());     \
  }

#define SET_PARAMETER_HELPER_VECTOR2D(PARAMETER_NAME, MEMBER_NAME)             \
  if (name == PARAMETER_NAME) {                                                \
    MEMBER_NAME = Vector2d(boost::get<std::vector<double>>(value).data());     \
  }

/**
 * Convinience typedefs.
 */
typedef std::map<std::string, Variant> VariantMap;
typedef std::map<std::string, Parameter> ParameterMap;

/**
 * @brief Make a Variant from argument.
 *
 * This is a convinience function, so that
 * rather involved constructors from
 * boost::variant are not needed in the
 * script interfaces.
 */
template <typename T> Variant make_variant(const T &x) { return Variant(x); }

/**
 * @brief Base class for generic script interface.
 *
 * @TODO Add extensive documentation.
 *
 */
class ScriptInterfaceBase : public Utils::ObjectId<ScriptInterfaceBase> {
public:
  /**
   * @brief Name of the object.
   *
   * Should be the name of the derived type including
   * namespace qualifiers. Must be unique.
   *
   * @return Name of the object.
   */
  // return boost::core::demangle(typeid(*this).name()) ?
  virtual const std::string name() const = 0;

  /**
   * @brief get current parameters.
   * @return Parameters set in class.
   */
  virtual VariantMap get_parameters() const {
    VariantMap values;

    for (auto const &p : valid_parameters()) {
      values[p.first] = get_parameter(p.first);
    }

    return values;
  }

  /**
   * @brief Get requiered and optional parameters for class
   *
   * Get requiered and optional parameters for class.
   *
   * @return Expected parameters.
   */
  virtual ParameterMap valid_parameters() const { return {}; }

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
  virtual void set_parameter(const std::string &name, const Variant &value){};
  /**
   * @brief Set multiple parameters.
   *
   * The default implementation calls the implementation of set_parameter for
   * every
   * element of the map.
   *
   * @param Paramters Parameters to set.
   */
  virtual void set_parameters(const VariantMap &parameters) {
    for (auto const &it : parameters) {
      set_parameter(it.first, it.second);
    }
  }

  /**
   * @brief Call a method on the object.
   *
   * If not overriden by the implementation,
   * this does nothing.
   */
  virtual Variant call_method(const std::string &, const VariantMap &) {}

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  static std::shared_ptr<ScriptInterfaceBase>
  make_shared(std::string const &name);

  /**
   * @brief Get a new reference counted instance of a script interface by
   * type.
   *
   */
  template <typename T> std::shared_ptr<T> static make_shared() {
    std::shared_ptr<T> sp = std::make_shared<T>();

    /* Id of the newly created instance */
    const int id = sp->id();

    /* Now get a reference to the corresponding weak_ptr in ObjectId and
       update
       it with our shared ptr, so that everybody uses the same ref count.
    */
    sp->get_instance(id) = std::static_pointer_cast<ScriptInterfaceBase>(sp);

    return sp;
  }
};

} /* namespace ScriptInterface */

#endif
