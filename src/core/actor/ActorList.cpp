/*
  Copyright (C) 2014-2018 The ESPResSo project

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

#include "ActorList.hpp"
#include <algorithm>
#include <cassert>

void ActorList::add(Actor *actor) { this->push_back(actor); }

void ActorList::remove(Actor *actor) {
  auto needle = std::find(this->begin(), this->end(), actor);
  assert(needle != this->end());
  this->erase(needle);
}
