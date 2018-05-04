from __future__ import print_function, absolute_import
include "myconfig.pxi"

IF SWIMMER_REACTIONS:
    cdef class Reaction:
        cdef _params
        cdef _ct_rate

    cdef extern from "swimmer_reaction.hpp":
        void reactions_sanity_checks()
        void local_setup_reaction()
        void integrate_reaction()

    cdef extern from "communication.hpp":
        void mpi_setup_reaction()
