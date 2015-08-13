#include "TclCommand.hpp"
#include <iostream>
#include <sstream>

#ifdef HAVE_CXX11
#include <type_traits>
#endif

using namespace std;

void TclCommand::add_subcommand(const TclCommand &c) {
  children.insert(pair<string, const TclCommand &>(c.m_impl.name(), c));
}

string TclCommand::print_to_string() {
  ostringstream res;

  res << m_impl.name();
  for(auto &p: m_impl.get_parameters()) {
    res << p.second.value;
  }

  return res.str();
}

void TclCommand::parse_from_string(list<string> &argv) {
  Parameters p = m_impl.get_parameters();

  for(list<string>::iterator it = argv.begin(); it != argv.end();) {
    string s = *it;
    cout << s << endl;
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
	for(int i = 0; i < p.n_elements; ++i) {
	  v.push_back(atof(it->c_str()));
	  it = argv.erase(it);
	}
	si->second = v;
	break;	
      }      
    }
  }
}

