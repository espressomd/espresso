#ifndef __TCLSCRIPTOBJECTMANAGER_HPP
#define __TCLSCRIPTOBJECTMANAGER_HPP

#include "TclCommand.hpp"
#include "ScriptObject.hpp"
#include "ObjectManager.hpp"
#include "Factory.hpp"

#ifdef HAVE_CXX11
#include <type_traits>
#endif

/** Tcl Interface for a ObjectManager<ScriptObject> */

template<class T>
class TclScriptObjectManager : public TclCommand {
public:
  TclScriptObjectManager(ObjectManager<T> &om, Tcl_Interp *_interp) : m_om(om), TclCommand(_interp) {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<ScriptObject, T>::value, "Type has to be subclass ob ScriptObject");
#endif
  }

  void parse_from_string(std::list<std::string> &argv) {
    switch(argv.size()) {
    case 0:
      /* @TODO: print all */
      break;
    case 1:
      {
        if(argv.front() == "new") {
          int id;
          stringstream ss(argv.front());
          ss >> id;
          if(ss.fail()) {
            throw string("usage()");
          }
          Tcl_AppendResult(interp, print_one(id).c_str(), 0);
          break;
        } else {
          Factory<T>::Instance().make(argv.front());
        }
      }
    }
  }

  std::string print_to_string() {
    return string("This should print all object in the container, but it does not.");
  }

private:
  ObjectManager<T> &m_om;

  string print_one(int id) {
    return TclScriptObject(m_om[id], interp).print_to_string();
  }
};

#endif
