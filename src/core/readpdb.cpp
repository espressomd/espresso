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
static void  add_lj_interaction(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &types, std::vector<PdbLJInteraction> interactions, const double rel_cutoff, const int first_type) {
  for(std::vector<PdbLJInteraction>::const_iterator it = interactions.begin(); it != interactions.end(); ++it) {
    for(std::set<PdbParser::itp_atomtype>::const_iterator jt = types.begin(); jt != types.end(); ++jt) {

      const double epsilon_ij = sqrt(it->epsilon * jt->epsilon);
      const double sigma_ij = 0.5*(it->sigma+10.*jt->sigma);
      const double cutoff_ij = rel_cutoff*sigma_ij;
      const double shift_ij = -(pow(sigma_ij/cutoff_ij,12) - pow(sigma_ij/cutoff_ij,6));
      READPDB_TRACE(printf("adding lj interaction types %d %d eps %e sig %e cut %e shift %e\n", it->other_type, first_type + jt.id, epsilon_ij, sigma_ij,
			   cutoff_ij, shift_ij););
      if((epsilon_ij <= 0) || (sigma_ij <= 0)) {
	continue;
      }
      else
	lennard_jones_set_params(it->other_type, first_type + jt->id, epsilon_ij, sigma_ij,
			       cutoff_ij, shift_ij, 0.0, -1.0, 0.0);
    }
  }
}

/* Add Lennard-Jones interactions between particles added from pdb/itp file */
static void add_lj_internal(std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &types, const double rel_cutoff, const int first_type) {
  for(std::set<PdbParser::itp_atomtype>::const_iterator it = types.begin(); it != types.end(); ++it) {
    for(std::set<PdbParser::itp_atomtype>::const_iterator jt = types.begin(); jt != types.end(); ++jt) {
      if(it->id > jt->id)
	continue;
      const double epsilon_ij = sqrt(it->epsilon * jt->epsilon);
      const double sigma_ij = 0.5*(10.*it->sigma+10.*jt->sigma);
      const double cutoff_ij = rel_cutoff*sigma_ij;
      const double shift_ij = -pow(sigma_ij/cutoff_ij,12) - pow(sigma_ij/cutoff_ij,6);
      READPDB_TRACE(printf("adding internal lj interaction types %d %d eps %e sig %e cut %e shift %e sigma_i %e sigma_j %e\n", first_type + it->id, first_type + jt->id, epsilon_ij, sigma_ij,
			   cutoff_ij, shift_ij););
      if((epsilon_ij <= 0) || (sigma_ij <= 0)) {
	continue;
      }
      else
	lennard_jones_set_params(first_type + it->id, first_type + jt->id, epsilon_ij, sigma_ij,
			       cutoff_ij, shift_ij, 0.0, -1.0, 0.0);      
    }
  }
}
#endif /* LENNARD_JONES */

static int add_particles(PdbParser::PdbParser &parser, int first_id, int default_type, std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> &seen_types, int first_type = 0, bool fit = false) {
  double pos[3];
  int id = first_id;
  int stat;
  int type;
  double q;
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
	const PdbParser::itp_atomtype &itp_atomtype = parser.itp_atomtypes[entry->second.type];
	seen_types.insert(itp_atomtype);
	type = itp_atomtype.id;
	q = entry->second.charge;
	READPDB_TRACE(printf("pdb-id %d es-id %d itp-type-id %d q %f es-type %d", it->i, id, type, q, first_type + type));
	READPDB_TRACE(std::cout << " res " << entry->second.type << std::endl);
      } else {	
	type = default_type;
	q = 0;
      }
      set_particle_type(id, first_type + type);
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
				char *itp_file, int first_type, bool fit, bool lj_internal) {
  int n_part;
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  if(itp_file) {
    if(!parser.parse_itp_file(itp_file))
      return 0;
  }

  /* Unique set of types that actually have particles */
  std::set<PdbParser::itp_atomtype, PdbParser::itp_atomtype_compare> seen_types;
  
  n_part = add_particles(parser, first_id, type, seen_types, first_type,fit);

#ifdef LENNARD_JONES
  add_lj_interaction(seen_types, ljInteractions, lj_rel_cutoff, first_type);
  if(lj_internal)
    add_lj_internal(seen_types, lj_rel_cutoff, first_type);
#endif

  return n_part;
}


