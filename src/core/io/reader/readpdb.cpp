#include "readpdb.hpp"
#include "grid.hpp"
#include "lj.hpp"

#include <iostream>
#include <cmath>
#include <set>

#ifdef READPDB_DEBUG
#define READPDB_TRACE(A) A
#else
#define READPDB_TRACE(A)
#endif

#ifdef LENNARD_JONES
/* Add user requested Lennard-Jones interactions */
static void  add_lj_interaction(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &types, std::vector<PdbLJInteraction> interactions, const double rel_cutoff) {
  for(std::vector<PdbLJInteraction>::const_iterator it = interactions.begin(); it != interactions.end(); ++it) {
    for(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare>::const_iterator jt = types.begin(); jt != types.end(); ++jt) {
      const double epsilon_ij = sqrt(it->epsilon * jt->epsilon);
      const double sigma_ij = 0.5*(it->sigma+10.*jt->sigma);
      const double cutoff_ij = rel_cutoff*sigma_ij;
      const double shift_ij = -(pow(sigma_ij/cutoff_ij,12) - pow(sigma_ij/cutoff_ij,6));
      if((epsilon_ij <= 0) || (sigma_ij <= 0)) {
	continue;
      }
      else {
	READPDB_TRACE(printf("adding lj interaction types %d %d eps %e sig %e cut %e shift %e\n", 
			   it->other_type, jt->espresso_id, epsilon_ij, sigma_ij, cutoff_ij, shift_ij););
	lennard_jones_set_params(it->other_type, jt->espresso_id, epsilon_ij, sigma_ij,
			       cutoff_ij, shift_ij, 0.0, 0.0);
      }
    }
  }
}

/* Add Lennard-Jones interactions between particles added from pdb/itp file */
static void add_lj_internal(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &types, const double rel_cutoff, bool only_diagonal) {
  for(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare>::const_iterator it = types.begin(); it != types.end(); ++it) {
    for(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare>::const_iterator jt = types.begin(); jt != types.end(); ++jt) {
      if(it->espresso_id > jt->espresso_id)
	continue;
      if(only_diagonal && (it->espresso_id != jt->espresso_id))
	continue;
      const double epsilon_ij = sqrtf(it->epsilon * jt->epsilon);
      const double sigma_ij = 0.5*(10.*it->sigma+10.*jt->sigma);
      const double cutoff_ij = rel_cutoff*sigma_ij;
      const double shift_ij = -pow(sigma_ij/cutoff_ij,12) - pow(sigma_ij/cutoff_ij,6);
      if((epsilon_ij <= 0) || (sigma_ij <= 0)) {
	continue;
      }
      else {
	READPDB_TRACE(printf("adding internal lj interaction types %d %d eps %e sig %e cut %e shift %e epsilon_i %e\n", 
			   it->espresso_id, jt->espresso_id, epsilon_ij, sigma_ij, cutoff_ij, shift_ij, it->epsilon););
	lennard_jones_set_params(it->espresso_id, jt->espresso_id, epsilon_ij, sigma_ij,
			       cutoff_ij, shift_ij, 0.0, 0.0);
      }
    }
  }
}
#endif /* LENNARD_JONES */

static int add_particles(PdbParser::PdbParser &parser, int first_id, int default_type, std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &seen_types, int first_type = 0, bool fit = false) {
  double pos[3];
  int id = first_id;
  int stat;
  int type;
#ifdef ELECTROSTATICS
  double q;
#endif
  PdbParser::BoundingBox bb;
  bb.llx = bb.lly = bb.llz = 0.0;
  double bb_l[3] = { box_l[0], box_l[1], box_l[2] };

  bb = parser.calc_bounding_box();

  if(fit) {
    bb_l[0] = (bb.urx - bb.llx);
    bb_l[1] = (bb.ury - bb.lly);
    bb_l[2] = (bb.urz - bb.llz);

    READPDB_TRACE(printf("bb dimensions (%f %f %f)\n", bb_l[0], bb_l[1], bb_l[2]));
    READPDB_TRACE(printf("bb ll (%f %f %f)\n", bb.llx, bb.lly, bb.llz));
    READPDB_TRACE(printf("bb ur (%f %f %f)\n", bb.urx, bb.ury, bb.urz));
    
    for(int i = 0; i < 3; i++) {
      if(bb_l[i] > box_l[i]) {
	rescale_boxl(i, bb_l[i]);
      }
    }
  }

  int last_type = first_type;
  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin(); it != parser.pdb_atoms.end(); ++it) {
    pos[0] = (it->x - bb.llx);
    pos[1] = (it->y - bb.lly);
    pos[2] = (it->z - bb.llz);

    stat = place_particle(id, pos);

    const std::map<int, PdbParser::itp_atom>::const_iterator entry = parser.itp_atoms.find(it->i);
    switch(stat) {
    case ES_PART_OK:
      std::cerr << "Warning: position and type of particle " << id << " was overwriten by value from pdb file." << std::endl;
      /* Fall through */
    case ES_PART_CREATED:
      /* See if we have a type from itp file, otherwise set default type */      
      if(entry != parser.itp_atoms.end()) {
	PdbParser::itp_atomtype itp_atomtype = parser.itp_atomtypes[entry->second.type];
	/* See if we have seen that type before */
	std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare>::iterator type_iterator = seen_types.find(itp_atomtype);
	if(type_iterator == seen_types.end()) {
	  itp_atomtype.espresso_id=last_type++;
	  type_iterator = seen_types.insert(itp_atomtype).first;
	}	
	itp_atomtype = *type_iterator;
#ifdef ELECTROSTATICS
	q = entry->second.charge;
#endif
	type = itp_atomtype.espresso_id;
	READPDB_TRACE(printf("pdb-id %d es-id %d itp-type-id %d q %f es-type %d", it->i, id, itp_atomtype.id, q, itp_atomtype.espresso_id));
	READPDB_TRACE(std::cout << " res " << entry->second.type << std::endl);
      } else {	
	type = default_type;
#ifdef ELECTROSTATICS
	q = 0;
#endif
      }
      set_particle_type(id, type);
#ifdef ELECTROSTATICS
      set_particle_q(id, q);
#endif
      id++;
      break;
    case ES_PART_ERROR:
      std::cerr << "Warning: Illegal particle id " << id << std::endl;
      return id - first_id;
      break;
    }
  }
  return id - first_id;
}

int pdb_add_particles_from_file(char *pdb_file, int first_id, int type, std::vector<PdbLJInteraction> &ljInteractions, double lj_rel_cutoff,
				char *itp_file, int first_type, bool fit, bool lj_internal, bool lj_diagonal) {
  int n_part;
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  if(itp_file) {
    if(!parser.parse_itp_file(itp_file))
      return 0;
  }

#ifdef READPDB_DEBUG
  for(std::map<int, PdbParser::itp_atom>::const_iterator it = parser.itp_atoms.begin(); it != parser.itp_atoms.end(); ++it) {
    printf("itp_atom id %d name '%s' q %e\n", it->second.i, it->second.type.c_str(), it->second.charge);
  }
#endif

  /* Unique set of types that actually have particles */
  std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> seen_types;
  
  n_part = add_particles(parser, first_id, type, seen_types, first_type,fit);

#ifdef READPDB_DEBUG
  for(std::map<std::string, PdbParser::itp_atomtype>::const_iterator it = parser.itp_atomtypes.begin(); it != parser.itp_atomtypes.end(); ++it) {
    printf("itp_atomtype id %d es_id %d name '%s' epsilon %e sigma %e\n", it->second.id, it->second.espresso_id, it->first.c_str(), it->second.epsilon, it->second.sigma);
  }
  for(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare>::const_iterator it = seen_types.begin(); it != seen_types.end(); ++it) {
    printf("atomtype %d espresso_type %d epsilon %e sigma %e\n", it->id, it->espresso_id, it->epsilon, it->sigma);
  }
#endif

#ifdef LENNARD_JONES
  add_lj_interaction(seen_types, ljInteractions, lj_rel_cutoff);
  if(lj_internal || lj_diagonal)
    add_lj_internal(seen_types, lj_rel_cutoff, lj_diagonal);
#endif

  return n_part;
}


