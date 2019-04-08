/*
   Copyright (C) 2019-2019 The ESPResSo project

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
/** \file
 *  Constants and enumerators for LB.
 */

#ifndef LB_CONSTANTS_HPP
#define LB_CONSTANTS_HPP

/** @brief Parameter fields for lattice Boltzmann
 *
 *  Determine what actions have to take place upon change of the respective
 *  parameter.
 */
enum class LBParam {
  DENSITY,   /**< fluid density */
  VISCOSITY, /**< fluid kinematic viscosity */
  AGRID,     /**< grid constant for fluid lattice */
  EXTFORCE,  /**< external force density acting on the fluid */
  BULKVISC,  /**< fluid bulk viscosity */
  KT         /**< thermal energy */
};

#endif /* LB_CONSTANTS_HPP */
