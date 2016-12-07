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
#include <string>
#include <fstream>
#include "Observable.hpp"
#include "integrate.hpp" 

namespace Observables {

Observable::Observable()
  : last_value(n_values())
  , n()
  , last_update(0)
  , autoupdate(0)
  , autoupdate_dt(0)
  , m_ofile()
  , m_filename()
  , m_binary(false)
{}


int Observable::calculate() {
  // Clear last value
  last_value.assign(n_values(), 0.0);

  int temp=actual_calculate();
  last_update = sim_time;
  return temp;
}


int Observable::update() {
  int temp=actual_update();
  last_update = sim_time;
  return temp;
}


void Observable::set_filename(std::string const& filename, bool binary)
{
  if ( !filename.empty() )
  {
    m_filename = filename;
    m_binary   = binary;

    auto mode = std::ios_base::app;
    if ( m_binary )
      mode |= std::ios_base::binary;

    if ( m_ofile.is_open() )
    {
      m_ofile.close();
      m_ofile.clear(); // clear flags
    }

    m_ofile.open(m_filename, mode);
  }
};


bool Observable::writable() const
{
  return m_ofile.is_open();
}


void Observable::write() const
{
  if ( writable() )
  {
    do_write();
  }
}


void Observable::do_write() const
{
  m_ofile << sim_time;
  for (auto p : last_value)
    m_ofile << " " << p;
  m_ofile << std::endl;
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
