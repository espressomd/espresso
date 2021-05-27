/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef _GLOBAL_HPP
#define _GLOBAL_HPP

/** @file
 *  This file contains the codes for global variables.
 *  Add any global variable that should be handled
 *  automatically in global.cpp. This includes
 *  distribution to other nodes and
 *  read/user-defined access from the script interface.
 */

/** @brief Check if all the global fields are synchronized between the nodes. */
void check_global_consistency();

/** @brief Field Enumeration
 */
enum Fields {
  /** index of \ref nptiso_struct::piston (only used for the events) */
  FIELD_NPTISO_PISTON,
};

#endif
