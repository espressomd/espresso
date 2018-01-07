/*
  Copyright (C) 2015,2016 The ESPResSo project

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

#ifndef UTILS_OBJECT_ID_HPP
#define UTILS_OBJECT_ID_HPP

#include <boost/serialization/access.hpp>

#include "NumeratedContainer.hpp"

#include <memory>
#include <string>
#include <type_traits>

namespace Utils {

template <typename T> class AutoObjectId;

template <typename T> class ObjectId {
public:
  ObjectId() : m_id(-1) {}

  bool operator==(ObjectId const &rhs) const { return m_id == rhs.m_id; }
  bool operator!=(ObjectId const &rhs) const { return m_id != rhs.m_id; }
  bool operator<(ObjectId const &rhs) const { return m_id < rhs.m_id; }

  int id() const { return m_id; }

  std::string to_string() const { return std::to_string(m_id); }

private:
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, unsigned int /* version */) {
    ar &m_id;
  }

  friend class AutoObjectId<T>;
  explicit ObjectId(unsigned i) : m_id(i) {}

  int m_id;
};

/**
 * @brief Class to automatically maintain a registry of
 * objects. The objects are assigned an id, by which an
 * instance can be retrieved.
 */
template <typename T> class AutoObjectId {
public:
  /* Assign an id on construction */
  AutoObjectId() : m_id(reg().add(std::weak_ptr<T>())) {}

  /* Remove id on destruction */
  virtual ~AutoObjectId() { reg().remove(m_id.m_id); }

  /**
   * @brief Get indentifier for this instance.
   */
  ObjectId<T> id() const { return m_id; }
  /**
   * @brief get instance by id. If the id is none (ObjectId()),
   * a nullptr is returned. If the id is unknown, an out-of-range
   * exception is thrown.
   */
  static std::weak_ptr<T> &get_instance(ObjectId<T> id) {
    return reg()[id.m_id];
  }

private:
  ObjectId<T> m_id;
  static Utils::NumeratedContainer<std::weak_ptr<T>> &reg() {
    static Utils::NumeratedContainer<std::weak_ptr<T>> m_reg(
        {{ObjectId<T>().id(), std::weak_ptr<T>()}});

    return m_reg;
  }
};

} /* namespace Utils */

#endif
