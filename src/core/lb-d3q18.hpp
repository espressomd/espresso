/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file lb-d3q18.hpp
 * Header file for the lattice Boltzmann D3Q18 model.
 *
 * This header file contains the definition of the D3Q18 model.
 */

#ifndef D3Q18_H
#define D3Q18_H

#ifdef LB

/** Velocity sub-lattice of the D3Q18 model */
static double d3q18_lattice[18][3] = { {  1.,  0.,  0. },
          			       { -1.,  0.,  0. },
          		               {  0.,  1.,  0. }, 
          			       {  0., -1.,  0. },
          			       {  0.,  0.,  1. }, 
          			       {  0.,  0., -1. },
          			       {  1.,  1.,  0. }, 
          			       { -1., -1.,  0. },
          			       {  1., -1.,  0. },
          			       { -1.,  1.,  0. },
          			       {  1.,  0.,  1. },
          			       { -1.,  0., -1. },
          			       {  1.,  0., -1. },
          			       { -1.,  0.,  1. },
          			       {  0.,  1.,  1. },
          			       {  0., -1., -1. },
          			       {  0.,  1., -1. },
          			       {  0., -1.,  1. } } ;

/** Coefficients for pseudo-equilibrium distribution of the D3Q18 model */
static double d3q18_coefficients[18][4] = { { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./12.,  1./6., 1./4., -1./4. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. },
					    { 1./24., 1./12., 1./8.,     0. } };

/** Coefficients in the functional for the equilibrium distribution */
static double d3q18_w[18] = { 1., 1., 1., 1., 1., 1.,
			      1./2., 1./2., 1./2., 1./2., 
			      1./2., 1./2., 1./2., 1./2.,
			      1./2., 1./2., 1./2., 1./2. };

LB_Model d3q18_model = { 18, d3q18_lattice, d3q18_coefficients, d3q18_w, 1./2. };

#endif /* LB */

#endif /* D3Q18_H */

/*@}*/
