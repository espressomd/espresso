/*
  Copyright (C) 2014,2015,2016 The ESPResSo project
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/* vim: set ts=8 sts=2 sw=2 et: */

#include "PdbParser.hpp"

#include <sstream>
#include <fstream>

using namespace std;

namespace PdbParser {

  BoundingBox PdbParser::calc_bounding_box() const {
    BoundingBox bb;

    bb.llx = std::numeric_limits<float>::max();
    bb.lly = std::numeric_limits<float>::max();
    bb.llz = std::numeric_limits<float>::max();
    bb.urx = -std::numeric_limits<float>::max();
    bb.ury = -std::numeric_limits<float>::max();
    bb.urz = -std::numeric_limits<float>::max();

    for(std::vector<pdb_atom>::const_iterator it = pdb_atoms.begin(); it != pdb_atoms.end(); ++it) {
      bb.llx = std::min(it->x, bb.llx);
      bb.lly = std::min(it->y, bb.lly);
      bb.llz = std::min(it->z, bb.llz);
      bb.urx = std::max(it->x, bb.urx);
      bb.ury = std::max(it->y, bb.ury);
      bb.urz = std::max(it->z, bb.urz);
    }  
    return bb;
  }

  bool PdbParser::parse_pdb_file(const string & filename) { 
    ifstream file;
    string tmp;
    pdb_atom a;
    
    pdb_atoms.clear();

    try {
      file.open(filename.c_str());
      while(file.good()) {

	file >> tmp;
	if(tmp == "ATOM") {
	  file.ignore(246,' ');
	  file >> a.i;
	  file >> tmp >> tmp >> tmp >> tmp;
	  file >> a.x >> a.y >> a.z;
	  pdb_atoms.push_back(a);
	}
      } 
    }
    catch (ifstream::failure& e) { // NOLINT
      return false;
    }

    return true; 
  }

  bool PdbParser::parse_itp_file(const string & filename) { 
    ifstream file(filename.c_str());
    string tmp, buf;
    itp_atom atom;
    std::size_t pos;    

    itp_atoms.clear();
    itp_atomtypes.clear();

    while(file.good()) {
      try {
	buf = char(file.get());
	/* Skipp leading whitespace */
	if(std::isspace(buf[0]))
	  continue;

	/* Comment, ignore rest of line */
	if(buf[0] == ';') {
	  std::getline(file, buf);
	  continue;
	}

	/* Section statement */
	if(buf == "[") {	  
	  std::getline(file, buf);
	  pos = buf.find_first_not_of(" \t[");
	  if(pos == std::string::npos)
	    continue;

	  std::string section = buf.substr(pos, std::string::npos);
	  pos = section.find_first_of(']');
	  section = section.substr(0, pos);
	  pos = section.find_last_not_of(" \t");
	  section = section.substr(0, pos+1);
	  
	  if(section == "atoms") {
	    while(file.good()) {
	      buf = char(file.get());

	      /* Ignore leading whitespace, check for end of file (standard says EOF is "generaly" -1) */
	      if(std::isspace(buf[0]) || (buf[0] == -1)) {
		continue;
	      }
	      /* End of atoms section */
	      if(buf[0] == '[') {
		file.unget();
		break;
	      }
	      /* Comment, ignore line */
	      if(buf[0] == ';') {
		std::getline(file, buf);
		continue;
	      }
	      /* Push back first char */
	      file.unget();
	      /* Parse line */
	      std::getline(file, buf);	      
	      std::istringstream line(buf);
	      line >> atom.i >> atom.type >> tmp >> tmp >> tmp >> tmp >> atom.charge;	      
	      itp_atoms.insert(std::pair<int, itp_atom>(atom.i, atom));
	    }
	  }
	  if(section == "atomtypes") {
	    itp_atomtype type;
	    std::string type_name;
	    while(file.good()) {
	      buf = char(file.get());

	      /* Ignore leading whitespace */
	      if(std::isspace(buf[0])) {
		continue;
	      }
	      /* End of atoms section */
	      if(buf[0] == '[') {
		file.unget();
		break;
	      }

	      /* Ignore leading whitespace, check for end of file (standard says EOF is "generaly" -1) */
	      if(std::isspace(buf[0]) || (buf[0] == -1)) {
		continue;
	      }

	      /* Push back first char */
	      file.unget();
	      /* Parse line */
	      std::getline(file, buf);
	      std::istringstream line(buf);
	      line >> type_name >> tmp >> tmp >> tmp >> tmp >> type.sigma >> type.epsilon;
	      /* Id is sequential number starting from zero */
	      type.id = itp_atomtypes.size();
	      itp_atomtypes.insert(std::pair<std::string, itp_atomtype>(type_name, type));
	    }
	  }
	}
      }
      catch (...) {
	return false;
      }
    }

    return true; 
  }

  bool PdbParser::parse_file(const string &  pdb_filename, const string &  itp_filename) {
    return parse_pdb_file(pdb_filename) && parse_itp_file(itp_filename);
  }

}
