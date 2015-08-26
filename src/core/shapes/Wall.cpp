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

#include "Wall.hpp"

using namespace std;

namespace Shapes {
  Parameters &Wall::all_parameters() const {
    static bool init = false;
    static Parameters p;
    if(!init) {
      p["normal"] = Parameter(Variant::DOUBLE_VECTOR, 3, true);
      p["d"] = Parameter(Variant::DOUBLE, true);
      init = true;
    }
    return p;
  }

  Parameters Wall::get_parameters() {
    Parameters p;

    p["normal"] = vector<double>(n, n+3);
    p["d"] = d;
      
    return p;
  }

  void Wall::set_parameters(Parameters &p) {
    for(Parameters::iterator it = p.begin(); it != p.end(); ++it) {
      if(!it->second.set)
        continue;

      if(it->first == "d")
        d = p["d"].value;

      if(it->first == "normal") {
        vector<double> v = p["normal"].value;        
        copy(v.begin(), v.end(), n);
      }
    }
  }

  int Wall::calculate_dist(const double *ppos, double *dist, double *vec)
  {
    int i;

    *dist = -d;
    for(i=0;i<3;i++) *dist += ppos[i]*n[i];
  
    for(i=0;i<3;i++) vec[i] = n[i] * *dist;
    return 0;
  }
}
