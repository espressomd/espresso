from .c_analyze cimport PartCfg, partCfg

cdef extern from "utils/checks/charge_neutrality.hpp" namespace "Utils":
    void check_charge_neutrality(PartCfg &partCfg)
