/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
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

#ifndef _RUNTIMEERRORCOLLECTOR_HPP
#define _RUNTIMEERRORCOLLECTOR_HPP

#include "config.hpp"
#include <list>
#include <string>
#include "mpi.h"

class RuntimeErrorCollector {
  std::list<std::string> errors;
  MPI_Comm comm;
public:
  RuntimeErrorCollector(MPI_Comm comm);

  void 
  warning(const std::string &msg,
          const char* function, const char* file, const int line);
  void 
  warning(const char *msg,
          const char* function, const char* file, const int line);
  void 
  warning(const std::ostringstream &mstr,
          const char* function, const char* file, const int line);

  void 
  error(const std::string &msg,
        const char* function, const char* file, const int line);
  void 
  error(const char *msg,
        const char* function, const char* file, const int line);
  void 
  error(const std::ostringstream &mstr,
        const char* function, const char* file, const int line);

  int count();
  void clear();

  std::list<std::string> gather();
  void gatherSlave();
};

#endif
