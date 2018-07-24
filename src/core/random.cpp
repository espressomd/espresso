/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010

    Max-Planck-Institute for Polymer Research, Theory Group

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

#include "random.hpp"
#include "communication.hpp"
#include "debug.hpp"

#include <sstream>

namespace Random {
using std::string;
using std::ostringstream;
using std::istringstream;
using std::vector;

using Communication::mpiCallbacks;

std::mt19937 generator;
std::normal_distribution<double> normal_distribution(0, 1);
std::uniform_real_distribution<double> uniform_real_distribution(0, 1);

bool user_has_seeded = false;

/** Local functions */

/**
 * @brief Get a string representation of the state of the PRNG.
 */
string get_state() {
  ostringstream os;
  os << generator;

  return os.str();
}

/**
 * @brief Set the state of the PRNG from a string representation.
 */
void set_state(const string &s) {
  istringstream is(s);
  is >> generator;
}

/**
 * @brief Get the state size of the PRNG
 */
int get_state_size_of_generator() {
  return generator.state_size; // this only works for the mersenne twister
                               // generator, other generators do not provide
                               // this member variable
}

/** Communication */

void mpi_random_seed_slave(int pnode, int cnt) {
  int this_seed;  user_has_seeded=true;

  MPI_Scatter(nullptr, 1, MPI_INT, &this_seed, 1, MPI_INT, 0, comm_cart);

  RANDOM_TRACE(printf("%d: Received seed %d\n", this_node, this_seed));
  init_random_seed(this_seed);
}

void mpi_random_seed(int cnt, vector<int> &seeds) {
  int this_seed;
  user_has_seeded = true;
  mpi_call(mpi_random_seed_slave, -1, cnt);

  MPI_Scatter(&seeds[0], 1, MPI_INT, &this_seed, 1, MPI_INT, 0, comm_cart);

  RANDOM_TRACE(printf("%d: Received seed %d\n", this_node, this_seed));

  init_random_seed(this_seed);
}

void mpi_random_set_stat_slave(int, int) {
  user_has_seeded = true;
  string msg;
  mpiCallbacks().comm().recv(0, SOME_TAG, msg);

  set_state(msg);
}

void mpi_random_set_stat(const vector<string> &stat) {
  user_has_seeded = true;
  mpi_call(mpi_random_set_stat_slave, 0, 0);

  for (int i = 1; i < n_nodes; i++) {
    mpiCallbacks().comm().send(i, SOME_TAG, stat[i]);
  }

  set_state(stat[0]);
}

void mpi_random_get_stat_slave(int, int) {
  string state = get_state();

  mpiCallbacks().comm().send(0, SOME_TAG, state);
}

string mpi_random_get_stat() {
  string res = get_state();

  mpi_call(mpi_random_get_stat_slave, 0, 0);

  for (int i = 1; i < n_nodes; i++) {
    string tmp;
    mpiCallbacks().comm().recv(i, SOME_TAG, tmp);
    res.append(" ");
    res.append(tmp);
  }

  return res;
}

void init_random(void) {
  /** Set the initial seed */
  init_random_seed(1 + this_node);

  /** Register callbacks */
  mpiCallbacks().add(mpi_random_seed_slave);
  mpiCallbacks().add(mpi_random_set_stat_slave);
  mpiCallbacks().add(mpi_random_get_stat_slave);
}

void init_random_seed(int seed)
{
  std::seed_seq seeder{seed}; //come up with "sane" initialization to avoid too many zeros in the internal state of the Mersenne twister
  generator.seed(seeder);
  generator.discard(1e6); //discard the first 1e6 random numbers to warm up the Mersenne-Twister PRNG
}

} /* Random */
