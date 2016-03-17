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

#ifndef __PARALLELOBJECT_HPP
#define __PARALLELOBJECT_HPP

#include <unordered_map>

namespace Utils {
  template<class T> class ParallelFactory;
}

class ParallelObject {
 public:
  ~ParallelObject();
  /**
   * @brief Get the global id of the instance.
   *
   * @return Global id of this object.
   */
  int get_id() const;

  /** This probably can and should be a pointer to const */
  static ParallelObject *get_local_address(int id);  
 protected:
  /** Run callback on slaves */
  void call_slaves(int par1, int par2);
  /**
   * @brief Default constructor.
   *
   * Construct a new ParallelObject instance.
   * This has to be run in parallel on all nodes
   * to insure proper registering of the callback and
   * local pointers. Object construction  should only
   * be run by the appropriat @class Utils::ParallelFactory,
   * which construacts an instance on all nodes and runs
   * the constructures in a clean way.
   */
  ParallelObject(); 
 private:
  /** Construction is only allowed via the ParallelFactory,
      otherwise there will be a MPI deadlock because the
      constructor tries to set up a callback on all nodes. */
  template<class T> friend class Utils::ParallelFactory;
  void *operator new(size_t size) {
    /** We don't want to change the functionality so we just call
        the default implementation */
    return ::operator new(size);
  }
  
  /** The callback function to run, if any. */
  virtual void callback(int par1, int par2) {};
  /** Id to identify the object and callback on all nodes */
  int m_callback_id;

  /** Mapping of global id to local instance */
  static std::unordered_map<int, ParallelObject *> address_table;
};

#endif
