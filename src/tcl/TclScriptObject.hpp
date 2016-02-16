#ifndef __TCLSCRIPTOBJECT_HPP
#define __TCLSCRIPTOBJECT_HPP

/** Tcl interface to ScriptObjects
 */

#include "TclCommand.hpp"
#include "ScriptObject.hpp"

class TclScriptObject : public TclCommand {
public:
  TclScriptObject(ScriptObject &so, Tcl_Interp* _interp) : m_so(so), TclCommand(_interp) {};
  std::string print_to_string();
  void parse_from_string(std::list<std::string> &argv);
  using TclCommand::create_command;
  void create_command() {
    create_command(m_so.name());
  }
  const std::string name() { return m_so.name(); };

private:
  ScriptObject &m_so;
};

#endif
