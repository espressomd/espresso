#ifndef __INTERACTION_CONSTRAINT_HPP
#define __INTERACTION_CONSTRAINT_HPP

#include "Constraint.hpp"
#include "particle_data.hpp"
#include "shapes/Shape.hpp"

#include <memory>

namespace Constraints {
enum class ReflectionType { NONE, NORMAL, NORMAL_TANGENTIAL };

class GeometryConstraint : public Constraint {
public:
  GeometryConstraint(int _penetrable = false,
                     ReflectionType _reflection_type = REFLECTION_NONE,
                     int _tuneable_slip = false, int _only_positive = 0,
                     int ia_type = -1)
      : penetrable(_penetrable), tuneable_slip(_tuneable_slip),
        reflection_type(_reflection_type) {
    part_rep.p.type = ia_type;
  }

  /** Public Methods */
  void add_energy(Particle *p, const double *folded_pos,
                  Observable_stat &energy);
  void add_force(Particle *p, const double *folded_pos);
  Shapes::Shape *get_shape() const { return m_shape; }
  virtual const std::string name() const {
    return Constraint::name() + std::string("GeometryConstraint::") +
           m_shape->name();
  }

  /** Data members */
  Particle part_rep;
  ReflectionType reflection_type;
  int penetrable;
  int tuneable_slip;
  int only_positive;

  /** Parsing stuff */
  std::map<std::string, Variant> get_parameters();
  Parameters all_parameters() const;
  void set_parameter(const std::string &name, const Variant &value);

private:
  /** Private methods */
  void reflect_particle(Particle *p, const double *distance_vector,
                        const double *folded_pos);

  /** Private data members */
  /** global id for shape, which is a ParallelObject */
  int m_shape_id;
  Shapes::Shape *m_shape;
};

} /* Constraints */

#endif
