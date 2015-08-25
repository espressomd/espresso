#ifndef __TCLCOMMAND_HPP
#define __TCLCOMMAND_HPP

#include "parser.hpp"

#include <list>
#include <string>

/** Base for something that can be registered with the Tcl interpreter.
 */

class TclCommand {  
public:
  TclCommand(Tcl_Interp* _interp) : interp(_interp) {};
  virtual void parse_from_string(std::list<std::string> &argv) = 0;
  virtual std::string print_to_string() = 0;
  void create_command(const std::string &command);
private:
  Tcl_Interp *interp;
};

#endif
