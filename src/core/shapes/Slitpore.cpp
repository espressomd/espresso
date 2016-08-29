/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#include "Slitpore.hpp"
/* For box_l */
#include "grid.hpp"

#include <cmath>

#define SQR(A) ((A)*(A))


using namespace std;

namespace Shapes {
  int Slitpore::calculate_dist(const double *ppos, double *dist, double *vec) {
    // the left circles
    double box_l_x = box_l[0];
    double c11[2] = { box_l_x/2-pore_width/2-upper_smoothing_radius, pore_mouth - upper_smoothing_radius };
    double c12[2] = { box_l_x/2-pore_width/2+lower_smoothing_radius, pore_mouth - pore_length  + lower_smoothing_radius };
    // the right circles
    double c21[2] = { box_l_x/2+pore_width/2+upper_smoothing_radius, pore_mouth - upper_smoothing_radius };
    double c22[2] = { box_l_x/2+pore_width/2-lower_smoothing_radius, pore_mouth - pore_length  + lower_smoothing_radius };

    //  printf("c11 %f %f\n", c11[0], c11[1]);
    //  printf("c12 %f %f\n", c12[0], c12[1]);
    //  printf("c21 %f %f\n", c21[0], c21[1]);
    //  printf("c22 %f %f\n", c22[0], c22[1]);


    if (ppos[2] > pore_mouth + channel_width/2) {
      //    printf("upper wall\n");
      // Feel the upper wall
      *dist = pore_mouth + channel_width - ppos[2];
      vec[0] = vec[1] = 0;
      vec[2] = -*dist;
      return 0;
    }

    if (ppos[0]<c11[0] || ppos[0] > c21[0]) {
      // Feel the lower wall of the channel
      //    printf("lower wall\n");
      *dist = ppos[2] - pore_mouth;
      vec[0] = vec[1] = 0;
      vec[2] = *dist;
      return 0;
    }

    if (ppos[2] > c11[1]) {
      // Feel the upper smoothing
      if (ppos[0] < box_l_x/2) {
        //    printf("upper smoothing left\n");
        *dist = sqrt( SQR(c11[0] - ppos[0]) + SQR(c11[1] - ppos[2])) - upper_smoothing_radius;
        vec[0] = -( c11[0] - ppos[0] ) * (*dist)/(*dist+upper_smoothing_radius);
        vec[1] = 0;
        vec[2] = -( c11[1] - ppos[2] ) * (*dist)/(*dist+upper_smoothing_radius);
        return 0;
      } else {
        //    printf("upper smoothing right\n");
        *dist = sqrt( SQR(c21[0] - ppos[0]) + SQR(c21[1] - ppos[2])) - upper_smoothing_radius;
        vec[0] = -( c21[0] - ppos[0] ) * (*dist)/(*dist+upper_smoothing_radius);
        vec[1] = 0;
        vec[2] = -( c21[1] - ppos[2] ) * (*dist)/(*dist+upper_smoothing_radius);
        return 0;
      }
    }
  
    if (ppos[2] > c12[1]) {
      // Feel the pore wall
      if (ppos[0] < box_l_x/2) {
        //    printf("pore left\n");
        *dist = ppos[0] - (box_l_x/2-pore_width/2);
        vec[0]=*dist;
        vec[1]=vec[2]=0;
        return 0;
      } else {
        //    printf("pore right\n");
        *dist =  (box_l_x/2+pore_width/2) - ppos[0];
        vec[0]=-*dist;
        vec[1]=vec[2]=0;
        return 0;
      }
    }

    if (ppos[0]>c12[0] && ppos[0] < c22[0]) {
      //    printf("pore end\n");
      // Feel the pore end wall
      *dist = ppos[2] - (pore_mouth-pore_length);
      vec[0]=vec[1]=0;
      vec[2]=*dist;
      return 0;
    }
    // Else
    // Feel the lower smoothing
    if (ppos[0] < box_l_x/2) {
      //    printf("lower smoothing left\n");
      *dist = -sqrt( SQR(c12[0] - ppos[0]) + SQR(c12[1] - ppos[2])) + lower_smoothing_radius;
      vec[0] = ( c12[0] - ppos[0] ) * (*dist)/(-*dist+lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c12[1] - ppos[2] ) * (*dist)/(-*dist+lower_smoothing_radius);
      return 0;
    } else {
      //    printf("lower smoothing right\n");
      *dist = -sqrt( SQR(c22[0] - ppos[0]) + SQR(c22[1] - ppos[2])) + lower_smoothing_radius;
      vec[0] = ( c22[0] - ppos[0] ) * (*dist)/(-*dist+lower_smoothing_radius);
      vec[1] = 0;
      vec[2] = ( c22[1] - ppos[2] ) * (*dist)/(-*dist+lower_smoothing_radius);
      return 0;
    }

    return 0;
  }
}

