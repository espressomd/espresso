from core.PartCfg cimport PartCfg

cdef extern from "statistics_chain.hpp":
    int chain_start
    int chain_n_chains
    int chain_length
    void calc_re(PartCfg &, double ** re)
    void calc_rg(PartCfg &, double ** rg)
    void calc_rh(PartCfg &, double ** rh)
