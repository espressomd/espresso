#include "GeometryConstraint.hpp"
#include "interaction_data.hpp"
#include "forces_inline.hpp"
#include "energy_inline.hpp"

using std::ostringstream;

namespace Constraints {
void GeometryConstraint::reflect_particle(Particle *p, const double *distance_vector, const double *folded_pos) {
  double vec[3];
  double norm; 

  memcpy(vec, distance_vector, 3*sizeof(double));

  norm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  p->r.p[0] = p->r.p[0]-2*vec[0];
  p->r.p[1] = p->r.p[1]-2*vec[1];
  p->r.p[2] = p->r.p[2]-2*vec[2];

  /* vec seams to be the vector that points from the wall to the particle*/
  /* now normalize it */ 
  switch(reflection_type) {
    case REFLECTION_NORMAL:
      vec[0] /= norm;
      vec[1] /= norm;
      vec[2] /= norm;
      /* calculating scalar product - reusing var norm */
      norm = vec[0] *  p->m.v[0] + vec[1] * p->m.v[1] + vec[2] * p->m.v[2];
      /* now add twice the normal component to the velcity */
      p->m.v[0] = p->m.v[0]-2*vec[0]*norm; /* norm is still the scalar product! */
      p->m.v[1] = p->m.v[1]-2*vec[1]*norm;
      p->m.v[2] = p->m.v[2]-2*vec[2]*norm;
      break;
    case REFLECTION_NORMAL_TANGENTIAL:
      /* if bounce back, invert velocity */
      p->m.v[0] =-p->m.v[0];
      p->m.v[1] =-p->m.v[1];
      p->m.v[2] =-p->m.v[2];
      break;  
    case REFLECTION_NONE:
      break;
  }
}

void GeometryConstraint::add_force(Particle *p, const double *folded_pos) {
  int j;
  double dist, vec[3], force[3], torque1[3], torque2[3];

  IA_parameters *ia_params;

  ia_params=get_ia_param(p->p.type, part_rep.p.type);
  dist=0.;
  for (j = 0; j < 3; j++) {
    force[j] = 0;
#ifdef ROTATION
    torque1[j] = torque2[j] = 0;
#endif
  }

  if(checkIfInteraction(ia_params)) {
    m_shape->calculate_dist(folded_pos, &dist, vec);

    if ( dist > 0 ) {
      calc_non_bonded_pair_force(p, &part_rep,
                                 ia_params,vec,dist,dist*dist, force,
                                 torque1, torque2);
#ifdef TUNABLE_SLIP
      if (tunable_slip) {
        add_tunable_slip_pair_force(p1, &constraints[n].part_rep,ia_params,vec,dist,force);
      }
#endif
    }
    else if (penetrable && (dist <= 0)) {
      if ( (!only_positive) && ( dist < 0 ) ) {
        calc_non_bonded_pair_force(p, &part_rep,
                                   ia_params,vec,-1.0*dist,dist*dist, force,
                                   torque1, torque2);
      }
    }
    else {
      if(reflection_type != REFLECTION_NONE){
        reflect_particle(p, vec, folded_pos);
      } else {
        ostringstream msg;
        msg <<"Constraint (" << name() << ")" << " violated by particle " << p->p.identity << " dist " << dist;
        runtimeError(msg);
      }
    }
  }
  for (j = 0; j < 3; j++) {
    p->f.f[j] += force[j];
    total_force[j] -= force[j];
#ifdef ROTATION
    p->f.torque[j] += torque1[j];
    part_rep.f.torque[j] += torque2[j];
#endif
  }
}

void GeometryConstraint::add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) {
  double dist, vec[3];
  IA_parameters *ia_params;
  double nonbonded_en = 0.0;

  ia_params = get_ia_param(p->p.type, part_rep.p.type);

  dist=0.;
  if(checkIfInteraction(ia_params)) {
    m_shape->calculate_dist(folded_pos, &dist, vec);
    if ( dist > 0 ) {
      nonbonded_en = calc_non_bonded_pair_energy(p, &part_rep, ia_params, vec, dist, dist*dist);
    }
    else if ( ( dist <= 0) && penetrable) {
      if ( !only_positive && (dist < 0) ) {
        nonbonded_en = calc_non_bonded_pair_energy(p, &part_rep,
                                                   ia_params, vec, -1.0*dist, dist*dist);
      }
    }
    else {
      ostringstream msg;
      msg <<"Constraint (" << name() << ")" << " violated by particle " << p->p.identity;
      runtimeError(msg);
    }
  }
  if (part_rep.p.type >= 0)
    *obsstat_nonbonded(&energy, p->p.type, part_rep.p.type) += nonbonded_en;
}

Parameters GeometryConstraint::get_parameters() {
  Parameters p = all_parameters();

  p["only_positive"] = only_positive;    
  p["type"] = part_rep.p.type;
  p["reflecting"] = reflection_type;
  p["penetrable"] = penetrable;
  p["tunable_slip"] = tuneable_slip;
  p["shape"] = m_shape_id;
    
  return p;
}

Parameters GeometryConstraint::all_parameters() const {
  Parameters p;

  p["only_positive"] = Parameter(Variant::INT, false);
  p["type"] = Parameter(Variant::INT, true);
  p["reflecting"] = Parameter(Variant::INT, false);
  p["penetrable"] = Parameter(Variant::INT, false);
  p["tunable_slip"] = Parameter(Variant::INT, false);
  p["shape"] = Parameter(Variant::INT, true);
      
  return p;
}

void GeometryConstraint::set_parameter(const std::string &name, const Variant &value) {
  std::cout << "GeometryConstraint::set_parameter(" << name << ", " << value << ")\n";
  SET_PARAMETER_HELPER("only_positive", only_positive);
  if(name == "type") {
    SET_PARAMETER_HELPER("type", part_rep.p.type);
    make_particle_type_exist(part_rep.p.type);
    return;
  }
  
  if(name == "reflecting")     
    reflection_type = static_cast<ReflectionType>(static_cast<int>(value));

  SET_PARAMETER_HELPER("penetrable", penetrable);
  SET_PARAMETER_HELPER("tunable_slip", tuneable_slip);

  if(name == "shape") {
    m_shape_id = value;
    m_shape = dynamic_cast<Shapes::Shape *>(ParallelObject::get_local_address(m_shape_id));
    assert(m_shape != nullptr);

    std::cout << "m_shape->name() = " << m_shape->name() << std::endl;
    
    return;
  }
}


}
