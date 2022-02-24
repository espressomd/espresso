/*
 * Copyright (C) 2014-2019 The ESPResSo project
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
#ifndef CORE_ACTOR_ACTOR_HPP
#define CORE_ACTOR_ACTOR_HPP

#include "SystemInterface.hpp"

/**
 * Generic abstract potential class.
 * All potentials should be derived from this one.
 */
class Actor {
public:
  virtual void computeForces(SystemInterface &) {}
  virtual void computeTorques(SystemInterface &) {}
  virtual void computeEnergy(SystemInterface &) {}
  virtual ~Actor() = default;
};

#endif /* CORE_ACTOR_ACTOR_HPP */
