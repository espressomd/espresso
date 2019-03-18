#ifndef ESPRESSO_OBSERVABLE_STAT_HPP
#define ESPRESSO_OBSERVABLE_STAT_HPP

#include <utils/List.hpp>

struct Observable_stat {
  /** Status flag for observable calculation.  For 'analyze energy': 0
      re-initialize observable struct, else everything is fine,
      calculation can start.  For 'analyze pressure' and 'analyze
      p_inst': 0 or !(1+v_comp) re-initialize, else all OK. */
  int init_status;

  /** Array for observables on each node. */
  Utils::List<double> data;

  /** number of Coulomb interactions */
  int n_coulomb;
  /** number of dipolar interactions */
  int n_dipolar;
  /** number of non bonded interactions */
  int n_non_bonded;
  /** Number of virtual sites relative (rigid body) contributions */
  int n_virtual_sites;
  /** Number of external field contributions */
  const static int n_external_field = 1;

  /** start of bonded interactions. Right after the special ones */
  double *bonded;
  /** start of observables for non-bonded interactions. */
  double *non_bonded;
  /** start of observables for Coulomb interaction. */
  double *coulomb;
  /** start of observables for Coulomb interaction. */
  double *dipolar;
  /** Start of observables for virtual sites relative (rigid bodies) */
  double *virtual_sites;
  /** Start of observables for external fields */
  double *external_fields;

  /** number of doubles per data item */
  int chunk_size;
};

/** Structure used only in the pressure and stress tensor calculation to
   distinguish
    non-bonded intra- and inter- molecular contributions. */
typedef struct {
  /** Status flag for observable calculation.
      For 'analyze energy': 0 re-initialize observable struct, else every thing
     is fine, calculation can start.
      For 'analyze pressure' and 'analyze p_inst': 0 or !(1+v_comp)
     re-initialize, else all OK. */
  int init_status_nb;

  /** Array for observables on each node. */
  Utils::List<double> data_nb;

  /** number of non bonded interactions */
  int n_nonbonded;

  /** start of observables for non-bonded intramolecular interactions. */
  double *non_bonded_intra;
  /** start of observables for non-bonded intermolecular interactions. */
  double *non_bonded_inter;

  /** number of doubles per data item */
  int chunk_size_nb;
} Observable_stat_non_bonded;

#endif // ESPRESSO_OBSERVABLE_STAT_HPP
