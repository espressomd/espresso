#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"
#include "shapes/NoWhere.hpp"
#include "shapes/Shape.hpp"

namespace Constraints {
class Constraint {
public:
  enum class ReflectionType { NONE, NORMAL, NORMAL_TANGENTIAL };

  Constraint()
      : m_shape(std::make_shared<Shapes::NoWhere>()),
        m_reflection_type(ReflectionType::NONE), m_penetrable(false),
        m_only_positive(false), m_tuneable_slip(0),
        m_ext_electric_field(0), m_ext_magn_field(0), m_type(-1){
    reset_force();
  }

  void add_energy(Particle *p, double *folded_pos,
                  Observable_stat &energy) const;

  void add_force(Particle *p, double *folded_pos);

  /* Calculate distance from the constraint */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }

  void set_shape(std::shared_ptr<Shapes::Shape> const &shape) {
    m_shape = shape;
  }

  Shapes::Shape const &shape() const { return *m_shape; }

  ReflectionType const &reflection_type() const;

  void reset_force() { m_local_force = Vector3d{0, 0, 0}; }
  int &only_positive() { return m_only_positive; }
  int &penetrable() { return m_penetrable; }
  double &ext_electric_field() { return m_ext_electric_field; }
  double &ext_magn_field() { return m_ext_magn_field; }
  int &type() { return m_type; }
  Vector3d total_force() const;

private:
  /** Private methods */
  void reflect_particle(Particle *p, const double *distance_vector,
                        const double *folded_pos) const;

  /** Private data members */
  std::shared_ptr<Shapes::Shape> m_shape;

  ReflectionType m_reflection_type;
  int m_penetrable;
  int m_only_positive;
  int m_tuneable_slip;
  double m_ext_electric_field;
  double m_ext_magn_field;
  int m_type;
  Vector3d m_local_force;
};

} /* namespace Constaints */

#endif
