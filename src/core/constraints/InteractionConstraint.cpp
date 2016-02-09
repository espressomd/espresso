#include "InteractionConstraint.hpp"
#include "forces_inline.hpp"
#include "energy_inline.hpp"

#include <iostream>

using std::ostringstream;

namespace Constraints {
  void InteractionConstraint::add_force(Particle *p, const double *folded_pos) {
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
      m_shape.calculate_dist(folded_pos, &dist, vec);

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

  void InteractionConstraint::add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) {
    std::cout << "InteractionConstraint::add_energy()" << std::endl;
    double dist, vec[3];
    IA_parameters *ia_params;
    double nonbonded_en = 0.0;

    ia_params = get_ia_param(p->p.type, part_rep.p.type);

    dist=0.;
    if(checkIfInteraction(ia_params)) {
      m_shape.calculate_dist(folded_pos, &dist, vec);
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

Parameters InteractionConstraint::get_parameters() {
  Parameters p = all_parameters();

  p["only_positive"] = only_positive;    
  p["type"] = part_rep.p.type;
  p["reflecting"] = reflection_type;
  p["penetrable"] = penetrable;
  p["tunable_slip"] = tuneable_slip;

  /** Add the parameters from the shape */
  Parameters shape_parameters = m_shape.get_parameters();
  for(auto &e: shape_parameters)
    p[e.first] = e.second;  
  
  for(auto &e: p)
    std::cout << e.first << " " << e.second.value << std::endl;
  
  return p;
}

Parameters InteractionConstraint::all_parameters() const {
  Parameters p;

  p["only_positive"] = Parameter(Variant::INT, false);
  p["type"] = Parameter(Variant::INT, true);
  p["reflecting"] = Parameter(Variant::INT, false);
  p["penetrable"] = Parameter(Variant::INT, false);
  p["tunable_slip"] = Parameter(Variant::INT, false);
    
  /** Add the parameters from the shape */
  Parameters shape_parameters = m_shape.all_parameters();
  p.insert(shape_parameters.begin(), shape_parameters.end());
    
  return p;
}

void InteractionConstraint::set_parameter(const std::string &name, const Variant &value) {
  std::cout << "InteractionConstraint::set_parameter(" << name << ", " << value << ")" << std::endl;
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
  
  m_shape.set_parameter(name, value);
}

}
