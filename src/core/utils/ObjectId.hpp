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

#include <memory>
#include <ostream>

#include <boost/serialization/access.hpp>

#include "NumeratedContainer.hpp"

namespace Utils {

template <typename T> class AutoObjectId {
public:
  class ObjectId {
  public:
    ObjectId() : m_id(-1) {}

    bool operator==(ObjectId const &rhs) const { return m_id == rhs.m_id; }
    bool operator!=(ObjectId const &rhs) const { return m_id != rhs.m_id; }
    bool operator<(ObjectId const &rhs) const { return m_id < rhs.m_id; }

  private:
    friend class AutoObjectId<T>;

    explicit ObjectId(unsigned i) : m_id(i) {}

    friend class boost::serialization::access;
    template <typename Archive>
    void serialize(Archive &ar, unsigned int /* version */) {
      ar &m_id;
    }

    int m_id;
  };

  /* Assign an id on construction */
  AutoObjectId() : m_id(reg().add(std::weak_ptr<T>())) {}

  /* Remove id on destruction */
  virtual ~AutoObjectId() { reg().remove(m_id.m_id); }

  /**
   * @brief Get indentifier for this instance.
   */
  ObjectId id() const { return m_id; }
  /**
   * @brief get instance by id.
   */
  static std::weak_ptr<T> &get_instance(ObjectId id) { return reg()[id.m_id]; }

private:
  ObjectId m_id;
  static Utils::NumeratedContainer<std::weak_ptr<T>> &reg() {
    static Utils::NumeratedContainer<std::weak_ptr<T>> m_reg;

    return m_reg;
  }
};

template <typename T> using ObjectId = typename AutoObjectId<T>::ObjectId;

} /* namespace Utils */

#endif
