#ifndef VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP

#include "config.hpp"
#ifdef VIRTUAL_SITES
#include "VirtualSites.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
/** @brief Virtual sites which are advected with an lb fuid. Forces on them are
 * instantaneously transferred to the fluid */
class VirtualSitesInertialessTracers : public VirtualSites {
  /** @brief Update positions and/or velocities of virtual sites

  * Velocities are only updated update_velocities() return true
  * @param recalc_positions can be used to skip the reculation of positions
  */
  void update(bool recalc_positions = true) const override{};
  /** Back-transfer forces (and torques) to non-virtual particles */
  void back_transfer_forces_and_torques() const override{};
  void after_force_calc() override;
  void after_lb_propagation() override;
  /** @brief Is a ghost communication needed after position updates */
  bool need_ghost_comm_after_pos_update() const override { return false; }
  /** Is a ghost comm needed before a velocity update */
  bool need_ghost_comm_before_vel_update() const override { return false; };
  /** Is a ghost comm needed before back_transfer */
  bool need_ghost_comm_before_back_transfer() const override { return false; };
};

#endif
#endif
#endif
