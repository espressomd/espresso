#ifndef LEES_EDWARDS_H
#define LEES_EDWARDS_H

#include "ParticleRange.hpp"
#include "algorithm/periodic_fold.hpp"
#include "boost/variant.hpp"
#include "config.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"

#include "utils/Vector.hpp"
/** \file lees_edwards.hpp
 *
 */

namespace LeesEdwards {

// Protocols determining shear rate and positional offset as a function of time
/** Lees Edwards protocol for un-altered periodic boundary conditions */
struct Off {
  Off() : m_shear_plane_normal_coord{0}, m_shear_dir_coord{0} {};
  double shear_velocity(double time) const { return 0.; };
  double pos_offset(double time) const { return 0.; };
  int m_shear_plane_normal_coord;
  int m_shear_dir_coord;
  static constexpr const bool verlet_list_support = true;
};

/** Lees-Edwards protocol for linear shearing */
struct LinearShear {
  LinearShear()
      : m_initial_pos_offset{0.}, m_shear_velocity{0.},
        m_shear_plane_normal_coord{0}, m_shear_dir_coord{0} {};
  double shear_velocity(double time) const { return m_shear_velocity; };
  double pos_offset(double time) const {
    return m_initial_pos_offset + time * m_shear_velocity;
  }
  double m_initial_pos_offset;
  double m_shear_velocity;
  int m_shear_plane_normal_coord;
  int m_shear_dir_coord;
  static constexpr const bool verlet_list_support = false;
};

/** Lees-Edwards protocol for osciallatory shearing */
struct OscillatoryShear {
  OscillatoryShear()
      : m_amplitude{0.}, m_frequency{0.}, m_time_0{0.},
        m_shear_plane_normal_coord{0}, m_shear_dir_coord{0} {};
  double pos_offset(double time) const {
    return m_amplitude * std::sin(m_frequency * (time - m_time_0));
  }
  double shear_velocity(double time) const {
    return m_frequency * m_amplitude *
           std::cos(m_frequency * (time - m_time_0));
  }
  double m_amplitude;
  double m_frequency;
  double m_time_0;
  int m_shear_plane_normal_coord;
  int m_shear_dir_coord;
  static constexpr const bool verlet_list_support = false;
};

/** Type which holds the currently active protocol */
using ActiveProtocol = boost::variant<Off, LinearShear, OscillatoryShear>;

/** Whether Lees-Edwards support the use of Verlet-Lists */
bool supports_verlet_list();

/** Visitor to get positional offset from the Lees-Edwards protocol */
class PosOffsetGetter : public boost::static_visitor<double> {
public:
  PosOffsetGetter(double time) : m_time{time} {}
  template <typename T> double operator()(const T &protocol) const {
    return protocol.pos_offset(m_time);
  }

private:
  double m_time;
};

inline double get_pos_offset(double time,
                             std::shared_ptr<ActiveProtocol> protocol) {
  if (!protocol)
    return 0.;
  return boost::apply_visitor(PosOffsetGetter(time), *protocol);
}

/** Visitor to get shear plane normal coord from the Lees-Edwards protocol */
class ShearPlaneNormalGetter : public boost::static_visitor<unsigned int> {
public:
  template <typename T> unsigned int operator()(const T &protocol) const {
    return protocol.m_shear_plane_normal_coord;
  }
};

inline unsigned int
get_shear_plane_normal_coord(std::shared_ptr<ActiveProtocol> protocol) {
  if (!protocol)
    return 0;
  return boost::apply_visitor(ShearPlaneNormalGetter(), *protocol);
}

/** Visitor to get shear direction coord from the Lees-Edwards protocol */
class ShearDirectionGetter : public boost::static_visitor<unsigned int> {
public:
  template <typename T> unsigned int operator()(const T &protocol) const {
    return protocol.m_shear_dir_coord;
  }
};

inline int get_shear_dir_coord(std::shared_ptr<ActiveProtocol> protocol) {
  if (!protocol)
    return 0;
  return boost::apply_visitor(ShearDirectionGetter(), *protocol);
}

/** Visitor to get shear velocity from the Lees-Edwards protocol */
class ShearVelocityGetter : public boost::static_visitor<double> {
public:
  ShearVelocityGetter(double time) : m_time{time} {};
  template <typename T> double operator()(const T &protocol) const {
    return protocol.shear_velocity(m_time);
  }

private:
  double m_time;
};

/** Calculation of current velocity*/
inline double get_shear_velocity(double time,
                                 std::shared_ptr<ActiveProtocol> protocol) {
  if (!protocol)
    return 0;
  return boost::apply_visitor(ShearVelocityGetter(time), *protocol);
}

/** Visitor to get Verlet list support from Lees-Edwards protocol */
class VerletListSupportGetter : public boost::static_visitor<bool> {
public:
  template <typename T> bool operator()(const T &protocol) const {
    return protocol.verlet_list_support;
  }

private:
  double m_time;
};

/** Calculation of current velocity*/
inline bool get_verlet_list_support(std::shared_ptr<ActiveProtocol> protocol) {
  if (!protocol)
    return true;
  return boost::apply_visitor(VerletListSupportGetter(), *protocol);
}

/** At the beginning of Lees_Edwards we have to reset all particle images
 * to zero*/
void local_image_reset(const ParticleRange &particles);

// Forward declaration
template <class BoxGeometry>
void push(Particle &p, double pos_offset, double shear_velocity,
          BoxGeometry &box) {
  if (!box.lees_edwards_protocol)
    return;
  auto const protocol = box.lees_edwards_protocol;
  auto const shear_dir = get_shear_dir_coord(protocol);
  auto const shear_plane_normal = get_shear_plane_normal_coord(protocol);
  int dummy_i;
  // Lees-Edwards acts at the zero step for the velocity half update and at
  // the midstep for the position.
  //
  // The update of the velocity at the end of the time step is triggered by
  // the flag and occurs in propagate_vel_finalize_p_inst
  if (p.r.p[shear_plane_normal] >= box.length()[shear_plane_normal]) {
    p.m.v[shear_dir] -= shear_velocity;
    p.r.p[shear_dir] -= pos_offset;
    p.p.lees_edwards_offset -= pos_offset;
    Algorithm::periodic_fold(p.r.p[shear_dir], &dummy_i,
                             box.length()[shear_dir]);
    p.p.lees_edwards_flag = -1; // perform a negative half velocity shift in
                                // propagate_vel_finalize_p_inst
  } else if (p.r.p[shear_plane_normal] < 0.) {
    p.m.v[shear_dir] += shear_velocity;
    p.r.p[shear_dir] += pos_offset;
    p.p.lees_edwards_offset += pos_offset;
    Algorithm::periodic_fold(p.r.p[shear_dir], &dummy_i,
                             box.length()[shear_dir]);
    p.p.lees_edwards_flag = 1; // perform a positive half velocity shift in
                               // propagate_vel_finalize_p_inst
  } else {
    p.p.lees_edwards_flag = 0;
  }
  p.p.lees_edwards_offset -=
      (0.5 * time_step * p.l.i[shear_plane_normal] * shear_velocity);
}

/** Type for cacheing the current Lees-Edwards position offset and shear
 * velocity */
struct Cache {
  Cache() : pos_offset{0.}, shear_velocity{0.} {};
  double pos_offset;
  double shear_velocity;
  void update(std::shared_ptr<ActiveProtocol> protocol, double pos_time,
              double vel_time) {
    pos_offset = LeesEdwards::get_pos_offset(pos_time, protocol);
    shear_velocity = LeesEdwards::get_shear_velocity(vel_time, protocol);
  }
};
} // namespace LeesEdwards
#endif
