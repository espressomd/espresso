#include "readpdb.hpp"

#include <iostream>

static int add_particles(PdbParser::PdbParser &parser, int first_id, int default_type, int first_type = 0) {
  double pos[3];
  int id = first_id;
  int stat;
  int type;
  for(std::vector<PdbParser::pdb_atom>::const_iterator it = parser.pdb_atoms.begin();
      it != parser.pdb_atoms.end(); ++it) {
    pos[0] = it->x;
    pos[1] = it->y;
    pos[2] = it->z;
    stat = place_particle(id, pos);

    const std::map<int, PdbParser::itp_atom>::const_iterator entry = parser.itp_atoms.find(it->i);

    switch(stat) {
    case ES_PART_OK:
      std::cerr << "Warning: position and type of particle " << id << " was overwriten by value from pdb file." << std::endl;
      /* Fall through */
    case ES_PART_CREATED:
      /* See if we have a type from itp file, otherwise set default type */      
      if(entry != parser.itp_atoms.end()) {
	type = parser.itp_atomtypes[entry->second.type].id;
      } else {	
	type = default_type;
      }
      set_particle_type(id, first_type + type);
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



int pdb_add_particles_from_file(char *pdb_file, int first_id, int type,
				char *itp_file, int first_type) {
  PdbParser::PdbParser parser;
  if(!parser.parse_pdb_file(pdb_file))
    return 0;

  if(itp_file) {
    if(!parser.parse_itp_file(itp_file))
      return 0;
  }
  
  return add_particles(parser, first_id, type, first_type);
}


