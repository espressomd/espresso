cdef extern from "communication.h":
    void mpi_bcast_parameter(int p)

cdef extern from "global.h":
    int FIELD_TEMPERATURE
    int FIELD_THERMO_SWITCH
    int FIELD_TEMPERATURE
    int FIELD_LANGEVIN_GAMMA

cdef extern from "thermostat.h":
    double temperature
    int thermo_switch
    double langevin_gamma
    int THERMO_OFF
    int THERMO_LANGEVIN
