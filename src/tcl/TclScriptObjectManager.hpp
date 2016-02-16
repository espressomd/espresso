#ifndef __TCLSCRIPTOBJECTMANAGER_HPP
#define __TCLSCRIPTOBJECTMANAGER_HPP

#include "TclCommand.hpp"
#include "ScriptObject.hpp"
#include "ObjectManager.hpp"

#ifdef HAVE_CXX11
#include <type_traits>
#endif

#include <iostream>
#include <string>
#include <exception>
#include <memory>

using std::string;

/** Tcl Interface for a ObjectManager<ScriptObject> */

template<class T>
class TclScriptObjectManager : public TclCommand {
public:
  TclScriptObjectManager(Tcl_Interp *_interp) : TclCommand(_interp) {
#ifdef HAVE_CXX11
    static_assert(std::is_base_of<ScriptObject, T>::value, "Type has to be subclass of ScriptObject");
#endif
  }

  void parse_from_string(std::list<std::string> &argv) {
    switch(argv.size()) {
    case 0:
      Tcl_AppendResult(interp, print_to_string().c_str(), 0);
      break;
    case 1:
      {
        int id;
        std::stringstream ss(argv.front());
        
        ss >> id;
        if(ss.fail()) {
          throw std::invalid_argument("Expected integer id");
        }        
        argv.pop_front();
        
        Tcl_AppendResult(interp, print_one(id).c_str(), 0);
        
        break;
      }
    default:
      {
        int id;
        if(argv.front() == "new") {
          argv.pop_front();
          id = m_om.add(argv.front());
          argv.pop_front();
        } else {
          std::stringstream ss(argv.front());
          ss >> id;
          if(ss.fail()) {
            throw std::invalid_argument("Expected integer id");
          }
          argv.pop_front();
        }
        
        TclScriptObject(*m_om[id], interp).parse_from_string(argv);
        std::stringstream ss;
        ss << id;
        Tcl_AppendResult(interp, ss.str().c_str(), 0);
        break;
      }
    }      
  }

  std::string print_to_string() {
    std::ostringstream ss;

    for(typename ObjectManager<T>::iterator it = m_om.begin(); it != m_om.end(); ++it) {
      ss << "{ " << print_one(it->first) << " } ";
    }

    return ss.str();    
  }


private:
  ObjectManager<T> m_om;

  string print_one(int id) {
    return std::string(m_om.name(id)).append(" ").append(TclScriptObject(*m_om[id], interp).print_to_string());
  }
};

#endif
