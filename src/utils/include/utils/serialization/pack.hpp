#ifndef UTILS_SERIALIZATION_PACK_HPP
#define UTILS_SERIALIZATION_PACK_HPP

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>

#include <sstream>
#include <string>

namespace Utils {
/**
 * @brief Pack a serialize type into a string.
 *
 * @tparam T Serializable type
 * @param v Value to serialize
 * @return String representation
 */
template <class T> std::string pack(T const &v) {
  std::stringstream ss;
  boost::archive::binary_oarchive(ss) << v;

  return ss.str();
}

/**
 * @brief Unpack a serialize type into a string.
 *
 * @tparam T Serializable type
 * @param state String to construct the value from, as returned by @function
 * pack.
 * @return Unpacked value
 */
template <class T> T unpack(std::string const &state) {
  namespace iostreams = boost::iostreams;

  iostreams::array_source src(state.data(), state.size());
  iostreams::stream<iostreams::array_source> ss(src);

  T val;
  boost::archive::binary_iarchive(ss) >> val;

  return val;
}
} // namespace Utils

#endif // UTILS_SERIALIZATION_PACK_HPP
