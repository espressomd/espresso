#ifndef LB_PARTICLE_COUPLING_HPP
#define LB_PARTICLE_COUPLING_HPP

#include "utils/Counter.hpp"

/** Calculate particle lattice interactions.
 *  So far, only viscous coupling with Stokesian friction is implemented.
 *  Include all particle-lattice forces in this function.
 *  The function is called from \ref force_calc.
 */
void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual);

uint64_t lb_lbcoupling_get_rng_state();
void lb_lbcoupling_set_rng_state(uint64_t counter);
void lb_lbcoupling_set_friction(double friction);
double lb_lbcoupling_get_friction();
void mpi_set_lb_coupling_counter_slave(int high, int low);

struct LB_Particle_Coupling {
  Utils::Counter<uint64_t> rng_counter_coupling;
  /*
   * @brief Friction constant for the particle coupling.
   */
  double friction = 0.0;

private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &rng_counter_coupling;
    ar &friction;
  }
};

extern LB_Particle_Coupling lb_particle_coupling;

#endif
