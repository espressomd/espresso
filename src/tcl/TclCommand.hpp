#ifndef __TCLCOMMAND_HPP
#define __TCLCOMMAND_HPP

#include "parser.hpp"

#include "ScriptObject.hpp"
#include <list>
#include <string>
#include <map>

class TclCommand {
public:
  TclCommand(ScriptObject &s, TclCommand *_p = 0) : m_impl(s) {};
  virtual void add_subcommand(const TclCommand &c);
  virtual void parse_from_string(std::list<std::string> &argv);
  virtual std::string print_to_string();
  
private:
  std::map<std::string, const TclCommand &> children;
  ScriptObject &m_impl;
};

/** Helper Function to get proper function pointer */
//template<class C>
//int TclScriptObject_wrapper(ClientData data, Tcl_Interp *interp, int argc, char *argv[]);

#endif
