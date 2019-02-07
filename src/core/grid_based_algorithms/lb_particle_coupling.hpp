#ifndef LB_PARTICLE_COUPLING_HPP
#define LB_PARTICLE_COUPLING_HPP

/** Calculate particle lattice interactions.
 *  So far, only viscous coupling with Stokesian friction is implemented.
 *  Include all particle-lattice forces in this function.
 *  The function is called from \ref force_calc.
 */
void calc_particle_lattice_ia();

uint64_t lb_coupling_rng_state();
void lb_coupling_set_rng_state(uint64_t counter);
void mpi_set_lb_coupling_counter_slave(int high, int low);

#endif
