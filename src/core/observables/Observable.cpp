/*
  Copyright (C) 2010,2011,2012,2013,2016 The ESPResSo project
  
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
#include "Observable.hpp"
#include "integrate.hpp" 

namespace Observables {

Observable::Observable() {
  last_update = 0;
  autoupdate = 0;
  autoupdate_dt = 0;
  last_value.resize(n_values());
};

int Observable::calculate() {
  int temp = 0;
  // Clear last value
  int n=n_values();
  last_value.resize(n);
  for (int i=0;i<last_value.size();i++)
    last_value[i]=0;

  temp=actual_calculate();
  last_update = sim_time;
  return temp;
}


int Observable::update() {
  int temp = 0;
  temp=actual_update();
  last_update = sim_time;
  return temp;
}



	


void autoupdate_observables() {
//  int i;
//  for (auto iter=observables.begin(); iter!=observables.end(); iter++) {
//    Observable* o=iter->second;
//    if (o->autoupdate && sim_time-o->last_update>o->autoupdate_dt*0.99999) {
//      o->update();
//    }
//  }
}

int observables_autoupdate=0;

} // Namespace Observables
