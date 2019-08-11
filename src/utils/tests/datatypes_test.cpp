/*
Copyright (C) 2019 The ESPResSo project

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
#define BOOST_TEST_MODULE abs test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Array.hpp>
#include <utils/Vector.hpp>
#include <utils/mpi/datatypes.hpp>
#include <utils/serialization/List.hpp>

struct NotAnMpiDatatype {};
struct AnMpiDatatype {};

namespace boost {
namespace mpi {
template <>
struct is_mpi_datatype<NotAnMpiDatatype> : public boost::mpl::false_ {};
template <> struct is_mpi_datatype<AnMpiDatatype> : public boost::mpl::true_ {};
} // namespace mpi
} // namespace boost

BOOST_AUTO_TEST_CASE(optional) {
  static_assert(
      not boost::mpi::is_mpi_datatype<boost::optional<NotAnMpiDatatype>>::value,
      "");
  static_assert(
      boost::mpi::is_mpi_datatype<boost::optional<AnMpiDatatype>>::value, "");
}

BOOST_AUTO_TEST_CASE(array) {
  static_assert(
      not boost::mpi::is_mpi_datatype<Utils::Array<NotAnMpiDatatype, 2>>::value,
      "");
  static_assert(
      boost::mpi::is_mpi_datatype<Utils::Array<AnMpiDatatype, 3>>::value, "");
}

BOOST_AUTO_TEST_CASE(vector) {
  static_assert(not boost::mpi::is_mpi_datatype<
                    Utils::Vector<NotAnMpiDatatype, 2>>::value,
                "");
  static_assert(
      boost::mpi::is_mpi_datatype<Utils::Vector<AnMpiDatatype, 3>>::value, "");
}

BOOST_AUTO_TEST_CASE(list) {
  static_assert(
      not boost::mpi::is_mpi_datatype<Utils::List<NotAnMpiDatatype>>::value,
      "");
  static_assert(boost::mpi::is_mpi_datatype<Utils::List<AnMpiDatatype>>::value,
                "");
}
