#ifndef VIRTUAL_SITES_VIRTUAL_SITES_OFF_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_OFF_HPP


#include "config.hpp" 
#ifdef VIRTUAL_SITES
#include "VirtualSites.hpp" 


/** @brief Do-nothing virtual-sites implementation */
   class VirtualSitesOff : public VirtualSites {
    /** @brief Update positions and/or velocities of virtual sites 

    * Velocities are only updated update_velocities() return true 
    * @param recalc_positions can be used to skip the reculation of positions 
    */
    void update(bool recalc_positions=true) const override {};
    /** Back-transfer forces (and torques) to non-virtual particles */
    void back_transfer_forces_and_torques() const override {};
    /** @brief Is a ghost communication needed after position updates */
    bool need_ghost_comm_after_pos_update() const override { return false;} 
    /** Is a ghost comm needed before a velocity update */
    bool need_ghost_comm_before_vel_update() const override {return false;};
    /** Is a ghost comm needed before back_transfer */
    bool need_ghost_comm_before_back_transfer() const override {return false;};
   };

#endif
#endif

