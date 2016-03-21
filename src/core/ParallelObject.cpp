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

#include <functional>

#include <boost/mpi.hpp>

#include "MpiCallbacks.hpp"
#include "ParallelObject.hpp"

std::unordered_map<int, ParallelObject *> ParallelObject::address_table;

/** Ctor and Dtor run on all nodes in parallel */
ParallelObject::ParallelObject() {
 {
    using namespace std::placeholders;    
    MpiCallbacks::function_type f;
    /** Bind member function to this instance */
    f = std::bind(&ParallelObject::callback, this, _1, _2);

    m_callback_id = mpiCallbacks().add(f);

    address_table[m_callback_id] = this;        
  }
}

ParallelObject::~ParallelObject() {
  /** Remove the callback when deleting the object */
  mpiCallbacks().remove(m_callback_id);
}

void ParallelObject::call_slaves(int par1, int par2) {
  mpiCallbacks().call(m_callback_id, par1, par2);
}
 
int ParallelObject::get_id() const { return m_callback_id; }

ParallelObject *ParallelObject::get_local_address(int id) {
  return address_table.at(id);
}
