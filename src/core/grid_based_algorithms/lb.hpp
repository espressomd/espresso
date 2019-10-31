/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
/** \file
 *
 *  %Lattice Boltzmann algorithm for hydrodynamic degrees of freedom.
 *
 *  For performance reasons it is clever to do streaming and collision at the
 *  same time because every fluid node has to be read and written only once.
 *  This increases mainly cache efficiency. This is achieved by
 *  @ref lb_collide_stream.
 *
 *  The hydrodynamic fields, corresponding to density, velocity and stress, are
 *  stored in @ref LB_FluidNode in the array @ref lbfields, the populations
 *  in @ref LB_Fluid in the array @ref lbfluid which is constructed as
 *  2 x (Nx x Ny x Nz) x 19 array.
 *
 *  Implementation in lb.cpp.
 */

#ifndef LB_H
#define LB_H

#include "config.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "grid_based_algorithms/lb-d3q19.hpp"
#include "grid_based_algorithms/lb_constants.hpp"

#include <array>
#include <boost/optional.hpp>
#include <memory>

#include "errorhandling.hpp"

#include "halo.hpp"

#include <utils/Counter.hpp>
#include <utils/Span.hpp>
#include <utils/serialization/multi_array.hpp>

/** Counter for the RNG */
extern boost::optional<Utils::Counter<uint64_t>> rng_counter_fluid;

/** Data structure for fluid on a local lattice site */
struct LB_FluidNode {
#ifdef LB_BOUNDARIES
  /** flag indicating whether this site belongs to a boundary */
  int boundary;
  Utils::Vector3d slip_velocity = {};
#endif // LB_BOUNDARIES

  /** local force density */
  Utils::Vector3d force_density;
#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
  // For particle update, we need the force on the nodes in LBM
  // Yet, Espresso resets the force immediately after the LBM update
  // Therefore we save it here
  Utils::Vector3d force_density_buf;
#endif
};

/** Data structure holding the parameters for the Lattice Boltzmann system. */
struct LB_Parameters {
  /** number density (LB units) */
  double density;

  /** kinematic viscosity (LB units) */
  double viscosity;

  /** bulk viscosity (LB units) */
  double bulk_viscosity;

  /** lattice spacing */
  double agrid;

  /** time step for fluid propagation (MD units)
   *  Note: Has to be larger than MD time step! */
  double tau;

  /** external force density applied to the fluid at each lattice site (LB
   * Units) */
  Utils::Vector3d ext_force_density;

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
  /*@{*/
  /** amplitudes of the fluctuations of the modes */
  Utils::Vector19d phi;
  /*@}*/
  /** Thermal energy */
  double kT;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &density &viscosity &bulk_viscosity &agrid &tau &ext_force_density
        &gamma_odd &gamma_even &gamma_shear &gamma_bulk &is_TRT &phi &kT;
  }
};

/** %Lattice Boltzmann parameters. */
extern LB_Parameters lbpar;

/** The underlying lattice */
extern Lattice lblattice;

/** Communicator for halo exchange between processors */
extern HaloCommunicator update_halo_comm;

void lb_init(const LB_Parameters &lb_parameters);

void lb_reinit_fluid(std::vector<LB_FluidNode> &lb_fields,
                     const Lattice &lb_lattice,
                     const LB_Parameters &lb_parameters);

void lb_reinit_parameters(LB_Parameters &lb_parameters);
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

void lb_sanity_checks(const LB_Parameters &lb_parameters);

/** Sets the equilibrium distributions.
 *  @param index Index of the local site
 *  @param density local fluid density
 *  @param momentum_density local fluid flux density
 *  @param stress local fluid stress
 */
void lb_set_population_from_density_momentum_density_stress(
    Lattice::index_t index, double density,
    Utils::Vector3d const &momentum_density, Utils::Vector6d const &stress);

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
#endif

double lb_calc_density(std::array<double, 19> const &modes,
                       const LB_Parameters &lb_parameters);
Utils::Vector3d lb_calc_momentum_density(std::array<double, 19> const &modes,
                                         Utils::Vector3d const &force_density);
Utils::Vector6d lb_calc_stress(std::array<double, 19> const &modes,
                               Utils::Vector3d const &force_density,
                               const LB_Parameters &lb_parameters);

/** Calculation of hydrodynamic modes.
 *
 *  @param index number of the node to calculate the modes for
 *  @retval Array containing the modes.
 */
std::array<double, 19> lb_calc_modes(Lattice::index_t index,
                                     const LB_Fluid &lb_fluid);

/**
 * @brief Get the populations as a function of density, flux density and stress.
 * @param density fluid density
 * @param momentum_density       fluid flux density
 * @param stress      fluid stress
 * @return 19 populations (including equilibrium density contribution).
 */
Utils::Vector19d lb_get_population_from_density_momentum_density_stress(
    double density, Utils::Vector3d const &momentum_density,
    Utils::Vector6d const &stress);

inline Utils::Vector19d lb_get_population(Lattice::index_t index) {
  Utils::Vector19d pop{};
  for (int i = 0; i < D3Q19::n_vel; ++i) {
    pop[i] = lbfluid[i][index] + D3Q19::coefficients[i][0] * lbpar.density;
  }
  return pop;
}

inline void lb_set_population(Lattice::index_t index,
                              const Utils::Vector19d &pop) {
  for (int i = 0; i < D3Q19::n_vel; ++i) {
    lbfluid[i][index] = pop[i] - D3Q19::coefficients[i][0] * lbpar.density;
  }
}

uint64_t lb_fluid_get_rng_state();
void lb_fluid_set_rng_state(uint64_t counter);
void lb_prepare_communication(HaloCommunicator &halo_comm,
                              const Lattice &lb_lattice);

#ifdef LB_BOUNDARIES
/** Bounce back boundary conditions.
 * The populations that have propagated into a boundary node
 * are bounced back to the node they came from. This results
 * in no slip boundary conditions.
 *
 * [cf. Ladd and Verberg, J. Stat. Phys. 104(5/6):1191-1251, 2001]
 */
void lb_bounce_back(LB_Fluid &lbfluid, const LB_Parameters &lb_parameters,
                    const std::vector<LB_FluidNode> &lb_fields);

#endif /* LB_BOUNDARIES */

void lb_calc_fluid_momentum(double *result, const LB_Parameters &lb_parameters,
                            const std::vector<LB_FluidNode> &lb_fields,
                            const Lattice &lb_lattice);
void lb_collect_boundary_forces(double *result);
void lb_initialize_fields(std::vector<LB_FluidNode> &fields,
                          LB_Parameters const &lb_parameters,
                          Lattice const &lb_lattice);
void lb_on_param_change(LBParam param);

/*@}*/

#endif /* _LB_H */
