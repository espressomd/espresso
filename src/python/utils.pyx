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

