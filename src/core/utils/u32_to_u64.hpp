#ifndef UTILS_U32_TO_U64_HPP
#define UTILS_U32_TO_U64_HPP

#include <cinttypes>
#include <utility>

namespace Utils {
constexpr inline uint64_t u32_to_u64(uint32_t high, uint64_t low) {
  return (static_cast<uint64_t>(high) << 32) | static_cast<uint64_t>(low);
}

constexpr inline std::pair<uint32_t, uint32_t> u64_to_u32(uint64_t in) {
  return {static_cast<uint32_t>(in >> 32), static_cast<uint32_t>(in)};
}

} // namespace Utils

#endif // ESPRESSO_U32_TO_U64_HPP
