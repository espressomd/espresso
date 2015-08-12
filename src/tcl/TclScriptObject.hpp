#ifndef __TCLSCRIPTOBJECT_HPP
#define __TCLSCRIPTOBJECT_HPP

#include "ScriptObject.hpp"
#include <list>
#include <string>
#include <map>

class TclScriptObject : public ScriptObject {
public:
  virtual std::string &command_name() = 0;
  virtual void add_subcommand(TclScriptObject *c);
  virtual void parse_from_string(std::list<std::string> &argv);
private:
  std::map<std::string, TclScriptObject *> children;
};

#endif
