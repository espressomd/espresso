cdef extern from "communication.hpp":
    void mpi_bcast_parameter(int p)

cdef extern from "global.hpp":
    int FIELD_TEMPERATURE
    int FIELD_THERMO_SWITCH
    int FIELD_TEMPERATURE
    int FIELD_LANGEVIN_GAMMA

cdef extern from "thermostat.hpp":
    double temperature
    int thermo_switch
    double langevin_gamma
    int THERMO_OFF
    int THERMO_LANGEVIN
