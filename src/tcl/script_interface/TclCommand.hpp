#ifndef __TCLCOMMAND_HPP
#define __TCLCOMMAND_HPP

#include "parser.hpp"

#include <list>
#include <string>

/**
 * @brief Base for something that can be registered with the Tcl interpreter.
 */
class TclCommand {  
public:
  TclCommand(Tcl_Interp* interp)
      : m_interp(interp)
  {};

  virtual ~TclCommand() {}
  
  /** Parse arguments. */
  virtual void parse_from_string(std::list<std::string> &argv) = 0;
  /** Print to string. */
  virtual std::string print_to_string() const = 0;
  
  /** Register command with the Tcl interpreter.
   *
   * @param command The command name in Tcl.
   */
  void create_command(const std::string &command);
  /** Our Tcl interpreter */
  Tcl_Interp *interp() {
    return m_interp;
  }

private:
  Tcl_Interp *m_interp;
};

#endif
