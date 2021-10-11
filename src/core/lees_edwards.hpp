#ifndef LEES_EDWARDS_HPP
#define LEES_EDWARDS_HPP

#include "BoxGeometry.hpp"
#include "lees_edwards_protocol.hpp"

#include "Particle.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <stdexcept>

extern std::weak_ptr<LeesEdwards::ActiveProtocol> lees_edwards_active_protocol;

namespace LeesEdwards {
inline void update_offset(Particle &p, const BoxGeometry &box,
                          double time_step) {
  auto const &le = box.clees_edwards_bc();
  p.l.lees_edwards_offset -=
      0.5 * time_step * p.l.i[le.shear_plane_normal] * le.shear_velocity;
}

inline void push(Particle &p, const BoxGeometry &box, double time_step) {
  if (box.type() != BoxType::LEES_EDWARDS)
    return;
  auto const &le = box.clees_edwards_bc();

  // Lees-Edwards acts at the zero step for the velocity half update and at
  // the midstep for the position.
  //
  // The update of the velocity at the end of the time step is triggered by
  // the flag and occurs in propagate_vel_finalize_p_inst
  if (p.r.p[le.shear_plane_normal] >= box.length()[le.shear_plane_normal]) {
    p.l.lees_edwards_flag = -1; // perform a negative half velocity shift in
                                // propagate_vel_finalize_p_inst
  } else if (p.r.p[le.shear_plane_normal] < 0) {
    p.l.lees_edwards_flag = 1; // perform a positive half velocity shift in
                               // propagate_vel_finalize_p_inst
  } else {
    p.l.lees_edwards_flag = 0;
  }

  p.m.v[le.shear_direction] += p.l.lees_edwards_flag * le.shear_velocity;
  p.r.p[le.shear_direction] += p.l.lees_edwards_flag * le.pos_offset;
  p.l.lees_edwards_offset -= p.l.lees_edwards_flag * le.pos_offset;
  // TODO: clarify whether the fold is needed
  //  p.r.p[le.shear_direction] = Algorithm::periodic_fold(
  //      p.r.p[le.shear_direction], box.length()[le.shear_direction]);
  update_offset(p, box, time_step);
}

inline double velocity_shift(short int le_flag, const BoxGeometry &box) {
  if (box.type() != BoxType::LEES_EDWARDS)
    return 0.;

  return le_flag * box.clees_edwards_bc().shear_velocity;
}

inline Utils::Vector3d shear_direction(const BoxGeometry &box) {
  auto const dir = box.clees_edwards_bc().shear_direction;
  if (dir == 0)
    return {1., 0., 0.};
  if (dir == 1)
    return {0., 1., 0.};
  if (dir == 2)
    return {0., 0., 1.};
  throw std::runtime_error("coordinate out of range");
}

extern std::shared_ptr<ActiveProtocol> active_protocol;
extern double pos_offset_at_last_resort;

inline void update_pos_offset(const ActiveProtocol &protocol, BoxGeometry &box,
                              double time) {
  if (box.type() == BoxType::LEES_EDWARDS) {
    box.lees_edwards_bc().pos_offset = get_pos_offset(time, protocol);
  }
}

inline void update_shear_velocity(const ActiveProtocol &protocol,
                                  BoxGeometry &box, double time) {
  if (box.type() == BoxType::LEES_EDWARDS) {
    box.lees_edwards_bc().shear_velocity = get_shear_velocity(time, protocol);
  }
}

inline void on_resort(const BoxGeometry &box) {
  pos_offset_at_last_resort = box.clees_edwards_bc().pos_offset;
}

inline Utils::Vector3d verlet_list_offset(const BoxGeometry &box) {
  return shear_direction(box) * std::fabs(box.clees_edwards_bc().pos_offset -
                                          pos_offset_at_last_resort);
}

inline void on_simtime_change(BoxGeometry &box, double time) {
  if (box.type() != BoxType::LEES_EDWARDS)
    return;
  if (!active_protocol)
    return;
  update_pos_offset(*active_protocol, box, time);
  update_shear_velocity(*active_protocol, box, time);
}
} // namespace LeesEdwards
#endif
