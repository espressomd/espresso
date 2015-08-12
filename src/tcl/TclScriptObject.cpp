#include "TclScriptObject.hpp"
#include <iostream>
#include <sstream>

void TclScriptObject::add_subcommand(TclScriptObject *c) {
  children.insert(std::pair<std::string, TclScriptObject *>(c->command_name(), c));
}

void TclScriptObject::parse_from_string(std::list<std::string> &argv) {
  Parameters p = get_parameters();

  for(std::list<std::string>::iterator it = argv.begin(); it != argv.end();) {
    std::string s = *it;
    std::cout << s << std::endl;
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
	std::stringstream ss(*it);
	int i;
	ss >> i;
	if(ss.fail()) {
	  std::ostringstream error;
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
	std::stringstream ss(*it);
	double d;
	ss >> d;
	if(ss.fail()) {
	  std::ostringstream error;
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
	std::vector<int> v;
	for(int j = 1; j <= p.n_elements; ++j) {
	  {
	    std::stringstream ss(*it);
	    int i;
	    ss >> i;
	    if(ss.fail()) {
	      std::ostringstream error;
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
	std::vector<double> v;
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
