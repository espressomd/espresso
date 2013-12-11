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

cdef checkTypeOrExcept(x,t,msg):
  """Checks that x is of type t, otherwise throws ValueError with the message msg.
     If x is an array/list/tuple, all elements are checked.
     """
  # Check whether x is an array/list/tuple or a single value
  if hasattr(x, "__getitem__"): 
    for i in range(len(x)):
      if not isinstance(x[i], t):
         raise ValueError(msg)
  else:
    if not isinstance(x, t):
       raise ValueError(msg)

