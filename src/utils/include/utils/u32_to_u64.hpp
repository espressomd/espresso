/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
