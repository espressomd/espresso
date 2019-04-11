/*
  Copyright (C) 2019-2019 The ESPResSo project

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

/** \file
 *  Unit tests for the PdbParser class.
 */

#define BOOST_TEST_MODULE PdbParser test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>

#include "PdbParser.hpp"

BOOST_AUTO_TEST_CASE(parser) {
  // create input files
  const std::string pdb_filepath("PdbParser.pdb");
  const std::string itp_filepath("PdbParser.itp");
  std::ofstream file;
  file.open(pdb_filepath);
  file << "ATOM      1 O1   LIG A   1       1.000   1.000   1.000  1.00\n";
  file << "ATOM      2 C1   LIG A   1       2.000  -0.500   1.000  1.00";
  file.close();
  file.open(itp_filepath);
  file << "[ atomtypes ]\n";
  file << ";name   bond_type   mass   charge  ptype   sigma       epsilon\n";
  file << " c        c        12.010  0.7700   A     1.000e-01   2.000e-01\n";
  file << "\n";
  file << "[ atoms ]\n";
  file << ";   nr     type  resnr residue  atom   cgnr    charge     mass\n";
  file << "     1        O      1    LIG     O1      1   -0.7700       16\n";
  file << "     2        c      1    LIG     C1      2    0.7700    12.01";
  file.close();

  PdbParser::PdbParser topology{};
  // check file parser
  const bool success = topology.parse_file(pdb_filepath, itp_filepath);
  BOOST_CHECK(success);
  // check bounding box
  const auto box = topology.calc_bounding_box();
  BOOST_CHECK(box.llx == 1.0);
  BOOST_CHECK(box.lly == -0.5);
  BOOST_CHECK(box.llz == 1.0);
  BOOST_CHECK(box.urx == 2.0);
  BOOST_CHECK(box.ury == 1.0);
  BOOST_CHECK(box.urz == 1.0);
  // check itp atomtypes
  const auto O_atomtype = topology.itp_atomtypes["O"];
  BOOST_CHECK(O_atomtype.id == 0);
  BOOST_CHECK(O_atomtype.sigma == 0.0f);
  BOOST_CHECK(O_atomtype.epsilon == 0.0f);
  const auto c_atomtype = topology.itp_atomtypes["c"];
  BOOST_CHECK(c_atomtype.id == 1);
  BOOST_CHECK(c_atomtype.sigma == 0.1f);
  BOOST_CHECK(c_atomtype.epsilon == 0.2f);
  // check itp atoms
  const auto O_itp_atom = topology.itp_atoms[1];
  BOOST_CHECK(O_itp_atom.i == 1);
  BOOST_CHECK(O_itp_atom.type == "O");
  BOOST_CHECK(O_itp_atom.charge == -0.77f);
  const auto c_itp_atom = topology.itp_atoms[0];
  BOOST_CHECK(c_itp_atom.i == 0);
  BOOST_CHECK(c_itp_atom.type.empty());
  BOOST_CHECK(c_itp_atom.charge == 0.0f);
  // check pdb atoms
  const auto O_pdbatom = topology.pdb_atoms[0];
  BOOST_CHECK(O_pdbatom.i == 1);
  BOOST_CHECK(O_pdbatom.x == 1.0f);
  BOOST_CHECK(O_pdbatom.y == 1.0f);
  BOOST_CHECK(O_pdbatom.z == 1.0f);
  const auto c_pdbatom = topology.pdb_atoms[1];
  BOOST_CHECK(c_pdbatom.i == 2);
  BOOST_CHECK(c_pdbatom.x == 2.0f);
  BOOST_CHECK(c_pdbatom.y == -0.5f);
  BOOST_CHECK(c_pdbatom.z == 1.0f);
}
