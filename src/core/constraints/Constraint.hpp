#ifndef CONSTRAINTS_CONSTRAINT_HPP
#define CONSTRAINTS_CONSTRAINT_HPP

#include <memory>

#include "energy.hpp"
#include "particle_data.hpp"
#include "shapes/Shape.hpp"

namespace Constraints {
class Constraint {
public:
  enum class ReflectionType { NONE, NORMAL, NORMAL_TANGENTIAL };

  explicit Constraint(std::shared_ptr<Shapes::Shape> const &shape)
      : m_shape(shape), m_reflection_type(ReflectionType::NONE),
        m_penetrable(false), m_only_positive(false), m_tuneable_slip(0),
        m_type(-1) {}

  void add_energy(Particle *p, double *folded_pos,
                  Observable_stat &energy) const;

  void add_force(Particle *p, double *folded_pos);

  /* Calculate distance from the constraint */
  int calc_dist(const double *pos, double *dist, double *vec) const {
    return m_shape->calculate_dist(pos, dist, vec);
  }

  Shapes::Shape const &shape() const { return *m_shape; }

  ReflectionType const &reflection_type() const;

  void reset_force() { m_total_force = Vector3d{0, 0, 0}; }

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
  int m_type;
  Vector3d m_total_force;
};

} /* namespace Constaints */

#endif
