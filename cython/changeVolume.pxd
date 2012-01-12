cdef extern from "config.h":
	pass

cdef extern from "grid.h":
	cdef void rescale_boxl(int dir, double d_new)

