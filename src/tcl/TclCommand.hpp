#ifndef __TCLCOMMAND_HPP
#define __TCLCOMMAND_HPP

#include "parser.hpp"

#include "ScriptObject.hpp"
#include <list>
#include <string>
#include <map>

class TclCommand {  
public:
  TclCommand(ScriptObject &so, Tcl_Interp* _interp) : m_so(so), interp(_interp) {};
  virtual void add_subcommand(const TclCommand &c);
  virtual void add_subcommand(const std::string &command, const TclCommand &c);
  virtual void parse_from_string(std::list<std::string> &argv);
  virtual std::string print_to_string();
  virtual void create_command(const std::string &command);
  virtual TclCommand &create_command() {
    create_command(m_so.name());
    return *this;
  }
private:
  static std::map<std::string, const TclCommand &> children;
  ScriptObject &m_so;
  Tcl_Interp *interp;
};

#endif
