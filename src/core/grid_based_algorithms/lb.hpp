/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file
 *
 * %Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 * Implementation in lb.cpp.
 */

#ifndef LB_H
#define LB_H

#include "config.hpp"

#include "lattice.hpp"

void mpi_set_lb_fluid_counter(int high, int low);

#ifdef LB

#include "errorhandling.hpp"

#include "halo.hpp"
#include "utils.hpp"

#include <utils/Span.hpp>

/** \name Parameter fields for Lattice Boltzmann
 * The numbers are referenced in \ref mpi_bcast_lb_params
 * to determine what actions have to take place upon change
 * of the respective parameter.
 */
/*@{*/
#define LBPAR_DENSITY 0   /**< fluid density */
#define LBPAR_VISCOSITY 1 /**< fluid kinematic viscosity */
#define LBPAR_AGRID 2     /**< grid constant for fluid lattice */
#define LBPAR_TAU 3       /**< time step for fluid propagation */
/** friction coefficient for viscous coupling between particles and fluid */
#define LBPAR_FRICTION 4
#define LBPAR_EXTFORCE 5 /**< external force density acting on the fluid */
#define LBPAR_BULKVISC 6 /**< fluid bulk viscosity */

/** Note these are used for binary logic so should be powers of 2 */
#define LB_COUPLE_NULL 1
#define LB_COUPLE_TWO_POINT 2
#define LB_COUPLE_THREE_POINT 4

/*@}*/
/** Some general remarks:
 *  This file implements the LB D3Q19 method to Espresso. The LB_Model
 *  construction is preserved for historical reasons and might be removed
 *  soon. It is constructed as a multi-relaxation time LB, thus all populations
 *  are converted to modes, then collision is performed and transfered back
 *  to population space, where the streaming is performed.
 *
 *  For performance reasons it is clever to do streaming and collision at the
 *  same time because every fluid node has to be read and written only once.
 *  This increases mainly cache efficiency. Two alternatives are implemented:
 *  stream_collide and collide_stream.
 *
 *  The hydrodynamic fields, corresponding to density, velocity and stress, are
 *  stored in LB_FluidNodes in the array lbfields, the populations in lbfluid
 *  which is constructed as 2 x (Nx x Ny x Nz) x 19 array.
 */

/** Description of the LB Model in terms of the unit vectors of the
 *  velocity sub-lattice and the corresponding coefficients
 *  of the pseudo-equilibrium distribution
 */
template <size_t N_vel = 19> struct LB_Model {
  /** number of velocities */
  static const constexpr int n_veloc = static_cast<int>(N_vel);

  /** unit vectors of the velocity sublattice */
  std::array<std::array<double, 3>, N_vel> c;

  /** coefficients in the pseudo-equilibrium distribution */
  std::array<std::array<double, 4>, N_vel> coeff;

  /** weights in the functional for the equilibrium distribution */
  std::array<double, N_vel> w;

  /** basis of moment space */
  std::array<std::array<int, N_vel>, N_vel> e_ki;

  /** normalization factors for the moment basis */
  std::array<double, N_vel> w_k;

  /** speed of sound squared */
  double c_sound_sq;

  /** transposed basis of moment space */
  std::array<std::array<int, N_vel>, N_vel> e_ki_transposed;
};

/** Data structure for fluid on a local lattice site */
struct LB_FluidNode {
#ifdef LB_BOUNDARIES
  /** flag indicating whether this site belongs to a boundary */
  int boundary;
  Vector3d slip_velocity = {};
#endif // LB_BOUNDARIES

  /** local force density */
  Vector3d force_density;
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  // For particle update, we need the force on the nodes in LBM
  // Yet, Espresso resets the force immediately after the LBM update
  // Therefore we save it here
  Vector3d force_density_buf;
#endif
};

/** Data structure holding the parameters for the Lattice Boltzmann system. */
struct LB_Parameters {
  /** number density (LB units) */
  double rho;

  /** kinematic viscosity (LB units) */
  double viscosity;

  /** bulk viscosity (LB units) */
  double bulk_viscosity;

  /** lattice spacing */
  double agrid;

  /** time step for fluid propagation (MD units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** friction coefficient for viscous coupling (LJ units) */
  double friction;

  /** external force density applied to the fluid at each lattice site (LB
   * Units) */
  Vector3d ext_force_density;

  /** relaxation of the odd kinetic modes */
  double gamma_odd;
  /** relaxation of the even kinetic modes */
  double gamma_even;
  /** relaxation rate of shear modes */
  double gamma_shear;
  /** relaxation rate of bulk modes */
  double gamma_bulk;

  /** Flag determining whether lbpar.gamma_shear, gamma_odd, and gamma_even are
   *  calculated from lbpar.gamma_shear in such a way to yield a TRT LB with
   *  minimized slip at bounce-back boundaries
   */
  bool is_TRT;

  /** \name Derived parameters */
  /** Flag indicating whether fluctuations are present. */
  int fluct;
  /** amplitudes of the fluctuations of the modes */
  Vector<19, double> phi;
};

/** The DnQm model to be used. */
extern LB_Model<> lbmodel;

/** %Lattice Boltzmann parameters. */
extern LB_Parameters lbpar;

/** The underlying lattice */
extern Lattice lblattice;

extern HaloCommunicator update_halo_comm;

void lb_realloc_fluid();

/** Pointer to the velocity populations of the fluid.
 *  lbfluid contains pre-collision populations, lbfluid_post
 *  contains post-collision populations
 */
using LB_Fluid = std::array<Utils::Span<double>, 19>;
extern LB_Fluid lbfluid;

class LB_Fluid_Ref {
public:
  LB_Fluid_Ref(std::size_t index, const LB_Fluid &lb_fluid)
      : m_index(index), m_lb_fluid(lb_fluid) {}
  template <std::size_t I> const auto &get() const {
    return m_lb_fluid[I][m_index];
  }

private:
  const std::size_t m_index;
  const LB_Fluid &m_lb_fluid;
};

namespace Utils {

template <std::size_t I> auto get(const LB_Fluid_Ref &lb_fluid) {
  return lb_fluid.get<I>();
}

} // namespace Utils

/** Pointer to the hydrodynamic fields of the fluid */
extern std::vector<LB_FluidNode> lbfields;

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Update the lattice Boltzmann system for one time step.
 *  This function performs the collision step and the streaming step.
 *  If external force densities are present, they are applied prior to the
 *  collisions. If boundaries are present, it also applies the boundary
 *  conditions.
 */
void lattice_boltzmann_update();

void lb_sanity_checks();

/** Calculates the equilibrium distributions.
    @param index Index of the local site
    @param rho local fluid density
    @param j local fluid speed
    @param pi local fluid pressure
*/
void lb_calc_n_from_rho_j_pi(const Lattice::index_t index, const double rho,
                             const std::array<double, 3> &j,
                             const std::array<double, 6> &pi);

/** Calculates the coupling of MD particles to the LB fluid.
 * This function  is called from \ref force_calc. The force is added
 * to the particle force and the corresponding momentum exchange is
 * applied to the fluid.
 * Note that this function changes the state of the fluid!
 */
void calc_particle_lattice_ia();

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
Vector3d lb_lbfluid_get_interpolated_force(const Vector3d &p);
#endif

void lb_calc_local_fields(Lattice::index_t index, double *rho, double *j,
                          double *pi);

/** Calculation of hydrodynamic modes.
 *
 *  @param index number of the node to calculate the modes for
 *  @param mode output pointer to a double[19]
 */
std::array<double, 19> lb_calc_modes(Lattice::index_t index);

/** Calculate the local fluid density.
 * The calculation is implemented explicitly for the special case of D3Q19.
 * @param index the local lattice site (Input).
 * @param rho   local fluid density
 */
inline void lb_calc_local_rho(Lattice::index_t index, double *rho) {
  // unit conversion: mass density
  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_rho in " << __FILE__
                      << __LINE__ << ": CPU LB not switched on.";
    *rho = 0;
    return;
  }

  double avg_rho = lbpar.rho;

  *rho = avg_rho + lbfluid[0][index] + lbfluid[1][index] + lbfluid[2][index] +
         lbfluid[3][index] + lbfluid[4][index] + lbfluid[5][index] +
         lbfluid[6][index] + lbfluid[7][index] + lbfluid[8][index] +
         lbfluid[9][index] + lbfluid[10][index] + lbfluid[11][index] +
         lbfluid[12][index] + lbfluid[13][index] + lbfluid[14][index] +
         lbfluid[15][index] + lbfluid[16][index] + lbfluid[17][index] +
         lbfluid[18][index];
}

/** Calculate the local fluid momentum.
 *  The calculation is implemented explicitly for the special case of D3Q19.
 *  @param[in]  index  Local lattice site
 *  @param[out] j      Local fluid speed
 */
inline void lb_calc_local_j(Lattice::index_t index, double *j) {
  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_calc_local_j in " << __FILE__ << __LINE__
                      << ": CPU LB not switched on.";
    j[0] = j[1] = j[2] = 0;
    return;
  }

  j[0] = lbfluid[1][index] - lbfluid[2][index] + lbfluid[7][index] -
         lbfluid[8][index] + lbfluid[9][index] - lbfluid[10][index] +
         lbfluid[11][index] - lbfluid[12][index] + lbfluid[13][index] -
         lbfluid[14][index];
  j[1] = lbfluid[3][index] - lbfluid[4][index] + lbfluid[7][index] -
         lbfluid[8][index] - lbfluid[9][index] + lbfluid[10][index] +
         lbfluid[15][index] - lbfluid[16][index] + lbfluid[17][index] -
         lbfluid[18][index];
  j[2] = lbfluid[5][index] - lbfluid[6][index] + lbfluid[11][index] -
         lbfluid[12][index] - lbfluid[13][index] + lbfluid[14][index] +
         lbfluid[15][index] - lbfluid[16][index] - lbfluid[17][index] +
         lbfluid[18][index];
}

/** Calculate the local fluid fields.
 * The calculation is implemented explicitly for the special case of D3Q19.
 *
 * @param index   Index of the local lattice site.
 * @param rho     local fluid density
 * @param j       local fluid speed
 * @param pi      local fluid pressure
 */
void lb_calc_local_fields(Lattice::index_t index, double *rho, double *j,
                          double *pi);

#ifdef LB_BOUNDARIES
inline void lb_local_fields_get_boundary_flag(Lattice::index_t index,
                                              int *boundary) {

  if (!(lattice_switch & LATTICE_LB)) {
    runtimeErrorMsg() << "Error in lb_local_fields_get_boundary_flag in "
                      << __FILE__ << __LINE__ << ": CPU LB not switched on.";
    *boundary = 0;
    return;
  }

  *boundary = lbfields[index].boundary;
}
#endif

inline void lb_get_populations(Lattice::index_t index, double *pop) {
  for (int i = 0; i < lbmodel.n_veloc; ++i) {
    pop[i] = lbfluid[i][index] + lbmodel.coeff[i][0] * lbpar.rho;
  }
}

inline void lb_set_populations(Lattice::index_t index,
                               const Vector<19, double> &pop) {
  for (int i = 0; i < lbmodel.n_veloc; ++i) {
    lbfluid[i][index] = pop[i] - lbmodel.coeff[i][0] * lbpar.rho;
  }
}

uint64_t lb_coupling_rng_state_cpu();
void lb_coupling_set_rng_state_cpu(uint64_t counter);
uint64_t lb_fluid_rng_state_cpu();
void lb_fluid_set_rng_state_cpu(uint64_t counter);
void lb_prepare_communication();
#endif

#endif /* _LB_H */
