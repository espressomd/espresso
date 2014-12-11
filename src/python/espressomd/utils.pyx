#
# Copyright (C) 2013,2014 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
cdef extern from "stdlib.h":
  void free(void* ptr)
  void* malloc(size_t size)
  void* realloc(void* ptr, size_t size)
            

cdef IntList* create_IntList_from_python_object(obj):
  cdef IntList* il
  il=<IntList*> malloc(sizeof(IntList))
  init_intlist(il)
  
  alloc_intlist(il, len(obj))
  for i in range(len(obj)):
    il.e[i] = obj[i]
    print il.e[i]

  return il

cdef checkTypeOrExcept(x,n,t,msg):
  """Checks that x is of type t and that n values are given, otherwise throws ValueError with the message msg.
     If x is an array/list/tuple, the type checking is done on the elements, and
     all elements are checked.
     Integers are accepted when a float was asked for.
     """
  # Check whether x is an array/list/tuple or a single value
  if n>1:
    if hasattr(x, "__getitem__"): 
      for i in range(len(x)):
        if not isinstance(x[i], t):
          if not (t==float and isinstance(x[i],int)):
             raise ValueError(msg + " -- Item "+str(i)+" was of type "+type(x[i]).__name__)
    else:
      # if n>1, but the user passed a single value, also throw exception
      raise ValueError(msg+" -- A single value was given but "+str(n)+" were expected.")
  else:
    # N=1 and a single value
    if not isinstance(x, t):
      if not (t==float and isinstance(x,int)):
        raise ValueError(msg+" -- Got an "+type(x).__name__)

