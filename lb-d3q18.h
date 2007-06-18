/* $Id$
 *
 * This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
 * It is therefore subject to the ESPResSo license agreement which you
 * accepted upon receiving the distribution and by which you are
 * legally bound while utilizing this file in any form or way.
 * There is NO WARRANTY, not even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * You should have received a copy of that license along with this
 * program; if not, refer to http://www.espresso.mpg.de/license.html
 * where its current version can be found, or write to
 * Max-Planck-Institute for Polymer Research, Theory Group, 
 * PO Box 3148, 55021 Mainz, Germany. 
 * Copyright (c) 2002-2007; all rights reserved unless otherwise stated.
 */

/** \file lb-d3q18.h
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
