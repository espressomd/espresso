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
#include "RuntimeErrorCollector.hpp"
#include <sstream>
#include "communication.hpp"

using namespace std;

RuntimeErrorCollector::
RuntimeErrorCollector(MPI_Comm _comm) {
  this->comm = _comm;
}

void RuntimeErrorCollector::
warning(const string &msg,
        const char* function, const char* file, const int line) {
  ostringstream ostr;
  ostr << "{ WARNING: ";
  ostr << msg;
  ostr << " in function " << function << " (" << file << ":" << line 
       << ") on node " << this_node;
  ostr << " } ";
  errors.push_back(ostr.str());
}

void RuntimeErrorCollector::
warning(const char *msg,
        const char* function, const char* file, const int line) {
  this->warning(string(msg), function, file, line);
}

void RuntimeErrorCollector::
warning(const ostringstream &mstr,
        const char* function, const char* file, const int line) {
  this->warning(mstr.str(), function, file, line);
}

void RuntimeErrorCollector::
error(const string &msg,
      const char* function, const char* file, const int line) {
  ostringstream ostr;
  ostr << "{ ERROR: ";
  ostr << msg;
  ostr << " in function " << function << " (" << file << ":" << line 
       << ") on node " << this_node;
  ostr << " } ";
  errors.push_back(ostr.str());
}

void RuntimeErrorCollector::
error(const char *msg,
      const char* function, const char* file, const int line) {
  this->error(string(msg), function, file, line);
}

void RuntimeErrorCollector::
error(const ostringstream &mstr,
      const char* function, const char* file, const int line) {
  this->error(mstr.str(), function, file, line);
}


int RuntimeErrorCollector::
count() {
  int numMessages = this->errors.size();
  MPI_Allreduce(MPI_IN_PLACE, &numMessages, 1, MPI_INT, MPI_SUM, this->comm);
  return numMessages;
}

void RuntimeErrorCollector::clear() {
  errors.clear();
}

list<string> RuntimeErrorCollector::
gather() {
  int numMessages = this->count();

  // If no processor encountered an error, return
  if (numMessages == 0) return list<string>();

  list<string> allerrors = this->errors;

  // subtract the number of messages on the master, as they are not to be sent
  numMessages -= this->errors.size();

  MPI_Status status;
  int count;
  for (int i = 0; i < numMessages; ++i) {
    // get the next message
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, this->comm, &status);
    MPI_Get_count(&status, MPI_CHAR, &count);

    char buffer[count];
    MPI_Recv(buffer, count, MPI_CHAR, 
             MPI_ANY_SOURCE, MPI_ANY_TAG, this->comm, MPI_STATUS_IGNORE);
    
    allerrors.push_back(string());
    string &s = allerrors.back();
    s.assign(buffer, count);
  }

  this->clear();

  return allerrors;
}

void RuntimeErrorCollector::
gatherSlave() {
  // If no processor encountered an error, return
  if (this->count() == 0) return;
  
  // send all messages
  for (list<string>::iterator it = this->errors.begin();
       it != errors.end(); ++it) {
    MPI_Send(const_cast<char*>(it->data()), it->length(), MPI_CHAR, 0, 42, this->comm);
  }

  // finally empty the list
  this->clear();
}

