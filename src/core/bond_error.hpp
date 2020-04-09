/*
 * Copyright (C) 2010-2020 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef ESPRESSO_BOND_ERROR_HPP
#define ESPRESSO_BOND_ERROR_HPP

#include <utils/Span.hpp>

#include <stdexcept>

void bond_broken_error(int id, Utils::Span<const int> partner_ids);

/**
 * Exception indicating that a particle id
 * could not be resolved.
 */
struct BondResolutionError : std::exception {};

/**
 * Exception indicating that a bond type
 * was unknown.
 */
struct BondUnknownTypeError : std::exception {
  explicit BondUnknownTypeError(int type) : type(type) {}

  int type;
};

/**
 * Exception indicating that a bond with an
 * unexpected number of partners was encountered.
 */
struct BondInvalidSizeError : std::exception {
  explicit BondInvalidSizeError(int size) : size(size) {}

  int size;
};

#endif // ESPRESSO_BOND_ERROR_HPP
