/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA

#include "ek_reactions.hpp"

#include <algorithm>

namespace EK {

EKReactions<walberla::EKReactionBase> ek_reactions;

void perform_reactions() {
  if (ek_reactions.empty()) {
    return;
  }

  std::for_each(ek_reactions.begin(), ek_reactions.end(),
                [](auto const &reaction) { reaction->perform_reaction(); });
}

} // namespace EK

#endif // WALBERLA
