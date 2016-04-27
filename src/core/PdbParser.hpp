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

#ifndef __PDBPARSER_HPP
#define __PDBPARSER_HPP

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <limits>

namespace PdbParser {

  struct BoundingBox {
    float llx, lly, llz;
    float urx, ury, urz;
  };

  typedef struct {
    int i; // index
    int m; // model index
    float x,y,z;
  } pdb_atom;

  typedef struct {
    int i;
    std::string type;
    float charge;
  } itp_atom;

  typedef struct {
    int id, espresso_id;
    float sigma,epsilon;
  } itp_atomtype;

  struct itp_atomtype_compare {
    bool operator() (const itp_atomtype &a, const itp_atomtype &b) { return a.id < b.id; }
  };

  class PdbParser {
  public:
    bool parse_pdb_file(std::string filename);
    bool parse_itp_file(std::string filename);
    bool parse_file(std::string pdb_filename, std::string itp_filename);
    BoundingBox calc_bounding_box() const;
    std::vector<pdb_atom> pdb_atoms;
    std::map<int, itp_atom> itp_atoms;
    std::map<std::string, itp_atomtype> itp_atomtypes;
  };
};

#endif
