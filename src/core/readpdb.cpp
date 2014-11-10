#include "readpdb.hpp"
#include "grid.hpp"

#include <iostream>
#include <limits>
#include <cmath>

struct BoundingBox {
  float llx, lly, llz;
  float urx, ury, urz;
};

#ifdef LENNARD_JONES
static void  add_lj_interaction(PdbParser::PdbParser &parser, std::vector<PdbLJInteraction> interactions, const double rel_cutoff) {
  for(std::vector<PdbLJInteraction>::const_iterator it = interactions.begin(); it != interactions.end(); ++it) {
    for(std::map<std::string, PdbParser::itp_atomtype> jt = parser.itp_atomtypes.begin(); jt != parser.itp_atomtypes.end(); ++jt) {
      if(it->other_type == jt->id)
	continue;
      const double epsilon_ij = sqrt(it->epsilon * jt->second.epsilon);
      const double sigma_ij = 0.5*(it->sigma+jt->second.epsilon);
      lennard_jones_set_params(it->other_type, jt->second.id, epsilon_ij, sigma_ij,
			       rel_cutoff*sigma_ij);
    }
  }
}
#endif

static BoundingBox calc_bounding_box(PdbParser::PdbParser &parser) {
  BoundingBox bb;

  bb.llx = std::numeric_limits<float>::max();
  bb.lly = std::numeric_limits<float>::max();
  bb.llz = std::numeric_limits<float>::max();
  bb.urx = std::numeric_limits<float>::min();
  bb.ury = std::numeric_limits<float>::min();
  bb.urz = std::numeric_limits<float>::min();

  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin();
      it != parser.pdb_atoms.end(); ++it) {
    bb.llx = (it->x < bb.llx) ? it->x : bb.llx;
    bb.lly = (it->y < bb.lly) ? it->y : bb.lly;
    bb.llz = (it->z < bb.llz) ? it->z : bb.llz;
    bb.urx = (it->x > bb.urx) ? it->x : bb.urx;
    bb.ury = (it->y > bb.ury) ? it->y : bb.ury;
    bb.urz = (it->z > bb.urz) ? it->z : bb.urz;
  }  
  return bb;
}

static int add_particles(PdbParser::PdbParser &parser, int first_id, int default_type, int first_type = 0, bool fit = false) {
  double pos[3];
  int id = first_id;
  int stat;
  int type;
  double q;
  BoundingBox bb;
  bb.llx = bb.lly = bb.llz = 0.0;
  double scalex = 1.0, scaley = 1.0, scalez = 1.0;

  if(fit) {
    bb = calc_bounding_box(parser);
    scalex = box_l[0] / (bb.urx - bb.llx);
    scaley = box_l[1] / (bb.ury - bb.lly);
    scalez = box_l[2] / (bb.urz - bb.llz);
    std::cout << "llx " << bb.llx << std::endl;
    std::cout << "lly " << bb.lly << std::endl;
    std::cout << "llz " << bb.llz << std::endl;
    std::cout << "urx " << bb.urx << std::endl;
    std::cout << "ury " << bb.ury << std::endl;
    std::cout << "urz " << bb.urz << std::endl;
    std::cout << "scalex " << scalex << std::endl;
    std::cout << "scaley " << scaley << std::endl;
    std::cout << "scalez " << scalez << std::endl;
  }

  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin();
      it != parser.pdb_atoms.end(); ++it) {
    pos[0] = scalex * (bb.llx + it->x);
    pos[1] = scaley * (bb.lly + it->y);
    pos[2] = scalez * (bb.llz + it->z);
    stat = place_particle(id, pos);

    const std::map<int, PdbParser::itp_atom>::const_iterator entry = parser.itp_atoms.find(it->i);

    switch(stat) {
    case ES_PART_OK:
      std::cerr << "Warning: position and type of particle " << id << " was overwriten by value from pdb file." << std::endl;
      /* Fall through */
    case ES_PART_CREATED:
      /* See if we have a type from itp file, otherwise set default type */      
      if(entry != parser.itp_atoms.end()) {
	const PdbParser::itp_atom itp_atom = parser.itp_atomtypes[entry->second.type];
	type = itp_atom.id;
	charge = itp_atom.charge;
      } else {	
	type = default_type;
      }
      set_particle_type(id, first_type + type);
      #ifdef ELECTROSTATICS
      set_particle_q(id, charge);
      #enfif
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

int pdb_add_particles_from_file(char *pdb_file, int first_id, int type, std::vector<PdbLJInteraction> &ljInteractions,
				char *itp_file, int first_type, bool fit) {
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  if(itp_file) {
    if(!parser.parse_itp_file(itp_file))
      return 0;
  }

  return add_particles(parser, first_id, type, first_type,fit);
}


