#include "TclScriptObject.hpp"
#include <iostream>

using namespace std;

string TclScriptObject::print_to_string() {
  ostringstream res;

  for(auto &p: m_so->get_parameters()) {
      res << p.first << " " << p.second.value << " ";
  }

  return res.str();
}


void TclScriptObject::parse_from_string(list<string> &argv) {
  cout << "TclScriptObject::parse_from_string()" << endl;
  Parameters p = m_so->get_parameters();
  for(list<string>::iterator it = argv.begin(); it != argv.end();) {
    string s = *it;
    Parameters::iterator si = p.find(s);
    if(si == p.end()) {
      ++it;
      continue;
    }
    it = argv.erase(it);
    Parameter p = si->second;
    switch(p.type) {
    case Variant::NONE:
      si->second.set = true;
      break;
    case Variant::INT:
      {
	stringstream ss(*it);
	int i;
	ss >> i;
	if(ss.fail()) {
	  ostringstream error;
	  error << s << " expects one integer argument, but got '" << *it << "'";
	  throw(error.str());
	}
	else
	  si->second.value = i;	
      }
      si->second.set = true;
      it = argv.erase(it);
      break;
    case Variant::DOUBLE:
      {
	stringstream ss(*it);
	double d;
	ss >> d;
	if(ss.fail()) {
	  ostringstream error;
	  error << s << " expects one float argument, but got '" << *it << "'";
	  throw(error.str());
	}
	else
	  si->second.value = d;	
      }
      si->second.set = true;
      it = argv.erase(it);
      break;
    case Variant::STRING:
      si->second.value = *it;
      si->second.set = true;
      it = argv.erase(it);
      break;
    case Variant::INT_VECTOR:
      {
	vector<int> v;
	for(int j = 1; j <= p.n_elements; ++j) {
	  {
	    stringstream ss(*it);
	    int i;
	    ss >> i;
	    if(ss.fail()) {
	      ostringstream error;
	      error << s << " expects " << p.n_elements << " integer arguments, but argument " << j << " was '" << *it << "'";
	      throw(error.str());
	    }
	    else
	      v.push_back(i);

	    it = argv.erase(it);
	  }
	}
	si->second = v;
	break;
      }
    case Variant::DOUBLE_VECTOR:
      {
	vector<double> v;
	for(int j = 1; j <= p.n_elements; ++j) {
	  {
	    stringstream ss(*it);
	    double i;
	    ss >> i;
	    if(ss.fail()) {
	      ostringstream error;
	      error << s << " expects " << p.n_elements << " float arguments, but argument " << j << " was '" << *it << "'";
	      throw(error.str());
	    }
	    else
	      v.push_back(i);

	    it = argv.erase(it);
	  }
	}
	si->second = v;
	break;
      }
    }
  }
  /** Drop unset */
  for(Parameters::iterator it = p.begin(); it != p.end(); ++it)
    if(!it->second.set)
      p.erase(it);

  m_so->set_parameters(p);
}

