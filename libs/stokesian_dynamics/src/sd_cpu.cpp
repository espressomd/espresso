#include <vector>

#if THRUST_VERSION < 100908
#define host_pinned_memory_resource universal_host_pinned_memory_resource
#endif

#include "sd.hpp"
#include "stokesian_dynamics/sd_cpu.hpp"

/** This executes the Stokesian Dynamics solver on the CPU. "_host" refers to
 *  data being stored on the host (in THRUST terminology, as opposed to 
 *  "device" storage), which is ordinary RAM.
 *
 *  \param x_host particle positions
 *  \param f_host particle forces and torques (total 6 entries per particle)
 *  \param a_host particle radii
 *  \param n_part number of particles
 *  \param eta dynamic viscosity of the ambient fluid
 *  \param sqrt_kT_Dt square root of kT/Dt, necessary for thermalization
 *  \param offset integer simulation time (steps) used as seed for the RNG
 *  \param seed global seed for the RNG used for the whole simulation
 *  \param flg certain bits set in this register correspond to certain features activated
 */
std::vector<double> sd_cpu(std::vector<double> const &x_host,
                           std::vector<double> const &f_host,
                           std::vector<double> const &a_host,
                           std::size_t n_part, double eta, double sqrt_kT_Dt,
                           std::size_t offset, std::size_t seed, int flg) {
  sd::solver<policy::host, double> viscous_force{eta, n_part};
  return viscous_force.calc_vel(x_host, f_host, a_host, sqrt_kT_Dt, offset, seed, flg);
}
