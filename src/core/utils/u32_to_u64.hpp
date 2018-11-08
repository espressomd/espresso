#ifndef UTILS_U32_TO_U64_HPP
#define UTILS_U32_TO_U64_HPP

#include <utility>
#include <cinttypes>

namespace Utils {
    constexpr uint64_t u32_to_u64(uint32_t high, uint64_t low) {
        return (static_cast<uint64_t>(high) << 32) | static_cast<uint64_t>(low);
    }

    constexpr std::pair <uint32_t, uint32_t> u64_to_u32(uint64_t in) {
        return {static_cast<uint32_t>(in >> 32), static_cast<uint32_t>(in)};
    }

}

#endif //ESPRESSO_U32_TO_U64_HPP
