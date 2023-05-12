/*
 * Copyright (C) 2021-2022 The ESPResSo project
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
#ifndef CORE_LEES_EDWARDS_LEES_EDWARDS_HPP
#define CORE_LEES_EDWARDS_LEES_EDWARDS_HPP

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "lees_edwards/protocols.hpp"

#include <cmath>
#include <memory>

#include <iostream>
namespace LeesEdwards {
class UpdateOffset {
protected:
  LeesEdwardsBC const &m_le;

public:
  UpdateOffset(BoxGeometry const &box) : m_le{box.lees_edwards_bc()} {}

  void operator()(Particle &p, double pos_prefactor = 1.0) const {
    // Disabled as long as we do not use a two step LE update
  }
};

class Push : public UpdateOffset {
  BoxGeometry const &m_box;

public:
  Push(BoxGeometry const &box) : UpdateOffset{box}, m_box{box} {}

  void operator()(Particle &p, double pos_prefactor = 1.0) const {
    // Lees-Edwards acts at the zero step for the velocity half update and at
    // the midstep for the position.
    //
    // The update of the velocity at the end of the time step is triggered by
    // the flag and occurs in @ref propagate_vel_finalize_p_inst
    if (p.pos()[m_le.shear_plane_normal] >=
        m_box.length()[m_le.shear_plane_normal]) {
      p.lees_edwards_flag() = -1; // perform a negative half velocity shift
    } else if (p.pos()[m_le.shear_plane_normal] < 0) {
      p.lees_edwards_flag() = 1; // perform a positive half velocity shift
    } else {
      p.lees_edwards_flag() = 0;
    }

    auto const dir = static_cast<double>(p.lees_edwards_flag());
    p.v()[m_le.shear_direction] += dir * m_le.shear_velocity;
    p.pos()[m_le.shear_direction] += pos_prefactor * dir * m_le.pos_offset;
    p.lees_edwards_offset() -= pos_prefactor * dir * m_le.pos_offset;
    fold_position(p.pos(), p.image_box(), m_box);
    //    UpdateOffset::operator()(p,pos_prefactor);
  }
};

inline double velocity_shift(short int le_flag, BoxGeometry const &box) {
  if (box.type() != BoxType::LEES_EDWARDS)
    return 0.;

  return le_flag * box.lees_edwards_bc().shear_velocity;
}

inline Utils::Vector3d shear_direction(BoxGeometry const &box) {
  return Utils::unit_vector<double>(box.lees_edwards_bc().shear_direction);
}

inline Utils::Vector3d verlet_list_offset(BoxGeometry const &box,
                                          double pos_offset_at_last_resort) {
  if (box.type() == BoxType::LEES_EDWARDS) {
    return shear_direction(box) * std::fabs(box.lees_edwards_bc().pos_offset -
                                            pos_offset_at_last_resort);
  }
  return {};
}

/** @brief Get currently active Lees-Edwards protocol. */
std::weak_ptr<ActiveProtocol> get_protocol();

/** @brief Set a new Lees-Edwards protocol. */
void set_protocol(std::shared_ptr<ActiveProtocol> new_protocol);

/** @brief Delete the currently active Lees-Edwards protocol. */
void unset_protocol();

} // namespace LeesEdwards
#endif
