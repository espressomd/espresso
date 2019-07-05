from core.statistics cimport Observable_stat, Observable_stat_non_bonded

cdef extern from "energy.hpp":
    cdef Observable_stat total_energy
    cdef Observable_stat_non_bonded total_energy_non_bonded
    cdef void master_energy_calc()
    cdef void init_energies(Observable_stat * stat)
    double calculate_current_potential_energy_of_system()
