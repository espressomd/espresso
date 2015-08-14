#include "TclCommand.hpp"
#include <iostream>
#include <sstream>

#ifdef HAVE_CXX11
#include <type_traits>
#endif

using namespace std;

static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]);

void TclCommand::add_subcommand(const TclCommand &c) {
  add_subcommand(c.m_so->name(), c);
}

void TclCommand::add_subcommand(const std::string &command, const TclCommand &c) {
  children.insert(pair<string, const TclCommand &>(command, c));
}

string TclCommand::print_to_string() {
  ostringstream res;

  res << m_so->name() << " ";
  for(auto &p: m_so->get_parameters()) {
    if(p.second.set)
      res << p.first << " " << p.second.value << " ";
  }

  return res.str();
}

void TclCommand::parse_from_string(list<string> &argv) {
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
  m_so->set_parameters(p);
}
  
  void TclCommand::create_command(const std::string &command) {
    printf("TclCommand::create_command() this = %p\n", this);
    Tcl_CreateCommand(interp, command.c_str(), (Tcl_CmdProc *)TclCommand_wrapper, reinterpret_cast<ClientData>(this), NULL);    
  }

static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  printf("TclCommand_wrapper(data = %p, interp = %p, argc = %d, argv = %p)\n",
	 data, interp, argc, argv);
  list<string> args(argv + 1, argv + argc);

  TclCommand *p = reinterpret_cast<TclCommand *>(data);

  try {
    p->parse_from_string(args);
    if(!args.empty()) {
      throw std::string("Unknown argument '").append(args.front()).append("'");
    }
  } catch(std::string &err) {
    Tcl_AppendResult(interp, err.c_str(), 0);
    return TCL_ERROR;
  }

  Tcl_AppendResult(interp, p->print_to_string(), 0);
  
  return TCL_OK;
}
