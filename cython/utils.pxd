
cdef extern from "utils.h":
  struct IntList:
    pass
  void init_intlist(IntList *il)
  void alloc_intlist(IntList *il, int size)
  void realloc_intlist(IntList *il, int size)

