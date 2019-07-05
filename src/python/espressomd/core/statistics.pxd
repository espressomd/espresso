from core.PartCfg cimport PartCfg
from core.utils cimport List
from core.utils cimport Vector3d, Vector3i
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string

cdef extern from "statistics.hpp":
    int n_part
    int n_part_conf
    int n_configs

    ctypedef struct Observable_stat:
        int init_status
        List[double] data
        int n_coulomb
        int n_dipolar
        int n_non_bonded
        int n_virtual_sites
        double * bonded
        double * non_bonded
        double * coulomb
        double * dipolar
        double * virtual_sites
        double * external_fields

    ctypedef struct Observable_stat_non_bonded:
        pass

    cdef void calc_structurefactor(PartCfg & , int * p_types, int n_types, int order, double ** sf)
    cdef vector[vector[double]] modify_stucturefactor(int order, double * sf)
    cdef double mindist(PartCfg & , const List[int] & set1, const List[int] & set2)
    cdef double min_distance2(Vector3d & pos1, Vector3d & pos2)
    cdef List[int] nbhood(PartCfg & , const Vector3d & pos, double r_catch, const Vector3i & planedims)
    cdef double distto(PartCfg & , const Vector3d & pos, int pid)
    cdef double * obsstat_bonded(Observable_stat * stat, int j)
    cdef double * obsstat_nonbonded(Observable_stat * stat, int i, int j)
    cdef double * obsstat_nonbonded_inter(Observable_stat_non_bonded * stat, int i, int j)
    cdef double * obsstat_nonbonded_intra(Observable_stat_non_bonded * stat, int i, int j)
    cdef vector[double] calc_linear_momentum(int include_particles, int include_lbfluid)
    cdef vector[double] centerofmass(PartCfg & , int part_type)
    cdef int calc_cylindrical_average(
        PartCfg & , vector[double] center, vector[double] direction,
        double length, double radius, int bins_axial, int bins_radial,
        vector[int] types, map[string, vector[vector[vector[double]]]] & distribution)

    void calc_rdf(PartCfg & , vector[int] p1_types, vector[int] p2_types,
                  double r_min, double r_max, int r_bins, vector[double] rdf)

    void calc_rdf_av(PartCfg & , vector[int] p1_types, vector[int] p2_types,
                     double r_min, double r_max, int r_bins, vector[double] rdf,
                     int n_conf)

    void angularmomentum(PartCfg & , int p_type, double * com)

    void momentofinertiamatrix(PartCfg & , int p_type, double * MofImatrix)

    void analyze_append(PartCfg & )

    void calc_part_distribution(
        PartCfg &, int * p1_types, int n_p1, int * p2_types, int n_p2,
        double r_min, double r_max, int r_bins, int log_flag, double * low,
        double * dist)
