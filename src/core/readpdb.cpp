#include "readpdb.hpp"
#include "grid.hpp"
#include "lj.hpp"

#include <iostream>
#include <cmath>

#ifdef READPDB_DEBUG
#define READPDB_TRACE(A) A
#else
#define READPDB_TRACE(A)
#endif

#ifdef LENNARD_JONES
static void  add_lj_interaction(PdbParser::PdbParser &parser, std::vector<PdbLJInteraction> interactions, const double rel_cutoff, const int first_type) {
  for(std::vector<PdbLJInteraction>::const_iterator it = interactions.begin(); it != interactions.end(); ++it) {
    for(std::map<std::string, PdbParser::itp_atomtype>::const_iterator jt = parser.itp_atomtypes.begin(); jt != parser.itp_atomtypes.end(); ++jt) {

      const double epsilon_ij = sqrt(it->epsilon * jt->second.epsilon);
      const double sigma_ij = 0.5*(it->sigma+jt->second.epsilon);
      const double cutoff_ij = rel_cutoff*sigma_ij;
      const double shift_ij = -pow(sigma_ij/cutoff_ij,12) - pow(sigma_ij/cutoff_ij,6);
      READPDB_TRACE(printf("adding lj interaction types %d %d eps %e sig %e cut %e shift %e\n", it->other_type, first_type + jt->second.id, epsilon_ij, sigma_ij,
			   cutoff_ij, shift_ij););
      lennard_jones_set_params(it->other_type, first_type + jt->second.id, epsilon_ij, sigma_ij,
			       cutoff_ij, shift_ij, 0.0, -1.0, 0.0);
    }
  }
}
#endif

static int add_particles(PdbParser::PdbParser &parser, int first_id, int default_type, int first_type = 0, bool fit = false) {
  double pos[3];
  int id = first_id;
  int stat;
  int type;
  double q;
  PdbParser::BoundingBox bb;
  bb.llx = bb.lly = bb.llz = 0.0;
  double scale = 1;

  if(fit) {
    bb = parser.calc_bounding_box();
    scale = std::min(box_l[0] / (bb.urx - bb.llx), std::min(box_l[1] / (bb.ury - bb.lly), box_l[2] / (bb.urz - bb.llz)));
    READPDB_TRACE(std::cout << "llx " << bb.llx << std::endl;);
    READPDB_TRACE(std::cout << "lly " << bb.lly << std::endl;);
    READPDB_TRACE(std::cout << "llz " << bb.llz << std::endl;);
    READPDB_TRACE(std::cout << "urx " << bb.urx << std::endl;);
    READPDB_TRACE(std::cout << "ury " << bb.ury << std::endl;);
    READPDB_TRACE(std::cout << "urz " << bb.urz << std::endl;);
    READPDB_TRACE(std::cout << "scale " << scale << std::endl;);
    READPDB_TRACE(std::cout << "bb.urx - bb.llx " << bb.urx - bb.llx << std::endl;);
    READPDB_TRACE(std::cout << "bb.ury - bb.lly " << bb.ury - bb.lly << std::endl;);
    READPDB_TRACE(std::cout << "bb.urz - bb.llz " << bb.urz - bb.llz << std::endl;);
  }

  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin(); it != parser.pdb_atoms.end(); ++it) {
    pos[0] = scale * (it->x - bb.llx);
    pos[1] = scale * (it->y - bb.lly);
    pos[2] = scale * (it->z - bb.llz);
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
	type = itp_atomtype.id;
	q = entry->second.charge;
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
				char *itp_file, int first_type, bool fit) {
  int n_part;
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  if(itp_file) {
    if(!parser.parse_itp_file(itp_file))
      return 0;
  }

  n_part = add_particles(parser, first_id, type, first_type,fit);

#ifdef LENNARD_JONES
  add_lj_interaction(parser, ljInteractions, lj_rel_cutoff, first_type);
#endif

  return n_part;
}


