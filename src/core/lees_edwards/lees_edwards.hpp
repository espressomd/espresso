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

namespace LeesEdwards {
class UpdateOffset {
protected:
  LeesEdwardsBC const &m_le;
  double const m_half_time_step;

public:
  UpdateOffset(BoxGeometry const &box, double time_step)
      : m_le{box.clees_edwards_bc()}, m_half_time_step{0.5 * time_step} {}

  void operator()(Particle &p) const {
    p.lees_edwards_offset() -= m_half_time_step *
                               p.image_box()[m_le.shear_plane_normal] *
                               m_le.shear_velocity;
  }
};

class Push : public UpdateOffset {
  Utils::Vector3d const &m_box_l;

public:
  Push(BoxGeometry const &box, double time_step)
      : UpdateOffset{box, time_step}, m_box_l{box.length()} {}

  void operator()(Particle &p) const {
    // Lees-Edwards acts at the zero step for the velocity half update and at
    // the midstep for the position.
    //
    // The update of the velocity at the end of the time step is triggered by
    // the flag and occurs in @ref propagate_vel_finalize_p_inst
    if (p.pos()[m_le.shear_plane_normal] >= m_box_l[m_le.shear_plane_normal]) {
      p.lees_edwards_flag() = -1; // perform a negative half velocity shift
    } else if (p.pos()[m_le.shear_plane_normal] < 0) {
      p.lees_edwards_flag() = 1; // perform a positive half velocity shift
    } else {
      p.lees_edwards_flag() = 0;
    }

    auto const dir = static_cast<double>(p.lees_edwards_flag());
    p.v()[m_le.shear_direction] += dir * m_le.shear_velocity;
    p.pos()[m_le.shear_direction] += dir * m_le.pos_offset;
    p.lees_edwards_offset() -= dir * m_le.pos_offset;
    // TODO: clarify whether the fold is needed
    //  p.pos()[m_le.shear_direction] = Algorithm::periodic_fold(
    //      p.pos()[m_le.shear_direction], m_box_l()[m_le.shear_direction]);
    UpdateOffset::operator()(p);
  }
};

inline double velocity_shift(short int le_flag, BoxGeometry const &box) {
  if (box.type() != BoxType::LEES_EDWARDS)
    return 0.;

  return le_flag * box.clees_edwards_bc().shear_velocity;
}

inline Utils::Vector3d shear_direction(BoxGeometry const &box) {
  return Utils::unit_vector<double>(box.clees_edwards_bc().shear_direction);
}

inline void update_pos_offset(ActiveProtocol const &protocol, BoxGeometry &box,
                              double time) {
  assert(box.type() == BoxType::LEES_EDWARDS);
  box.lees_edwards_bc().pos_offset = get_pos_offset(time, protocol);
}

inline void update_shear_velocity(ActiveProtocol const &protocol,
                                  BoxGeometry &box, double time) {
  assert(box.type() == BoxType::LEES_EDWARDS);
  box.lees_edwards_bc().shear_velocity = get_shear_velocity(time, protocol);
}

inline Utils::Vector3d verlet_list_offset(BoxGeometry const &box,
                                          double pos_offset_at_last_resort) {
  if (box.type() == BoxType::LEES_EDWARDS) {
    return shear_direction(box) * std::fabs(box.clees_edwards_bc().pos_offset -
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
