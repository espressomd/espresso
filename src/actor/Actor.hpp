/*
  Copyright (C) 2014 The ESPResSo project

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
#ifndef _ACTOR_ACTOR_HPP
#define _ACTOR_ACTOR_HPP

#include "SystemInterface.hpp"

/**
 * Generic abstract potential class.
 * All potentials should be derived from this one.
 */
class Actor {
public:
	virtual void computeForces(SystemInterface &s) = 0;
	virtual ~Actor() {}
};

#endif /* _ACTOR_ACTOR_HPP */
