/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

#include "timer_tcl.hpp"

#include <string>
#include <iostream>
#include <sstream>

#include "utils/Timer.hpp"

using namespace std;

using Utils::Timing::Timer;

namespace {
string format_stats(const Utils::Timing::Timer::Stats &stats) {
  ostringstream ss;

  ss << stats.avg() << " " << stats.sig() << " " << stats.var() << " " << stats.min() << " " << stats.max() << " " << stats.n();
  return ss.str();
}
}

vector<map<string, Utils::Timing::Timer::Stats> >
mpi_gather_timers();

int tclcommand_timer(ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  auto timers = mpi_gather_timers();
  char buf[8];
  
  for(int i = 0; i < timers.size(); ++i) {
    snprintf(buf, 8, "%d", i);

    for(auto &it : timers[i]) {      
      Tcl_AppendResult(interp, "{ ", buf, " ", it.first.c_str(), " ", format_stats(it.second).c_str(), " } ", nullptr);
    }
  }

  return TCL_OK;
}
