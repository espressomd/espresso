#include "readpdb.hpp"

#include <iostream>

static int add_particles(PdbParser::PdbParser &parser, int first_id, int type) {
  double pos[3];
  int id = first_id;
  int stat;
  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin();
      it != parser.pdb_atoms.end(); ++it) {
    pos[0] = it->x;
    pos[1] = it->y;
    pos[2] = it->z;
    stat = place_particle(id, pos);

    switch(stat) {
    case ES_PART_OK:
      std::cerr << "Warning: position and type of particle " << id << " was overwriten by value from pdb file." << std::endl;
      /* Fall through */
    case ES_PART_CREATED:
      set_particle_type(id, type);
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

int pdb_add_particles_from_file(char *pdb_file, int first_id, int type) {
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  return add_particles(parser, first_id, type);
}

