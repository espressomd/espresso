#include "TclCommand.hpp"
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

/**
 * @brief Wraper function for tcl commands.
 *
 * The Tcl c interface needs a function pointer
 * as callback. This function takes a pointer to
 * a TclCommand instance and calls the apopriate
 * methods on it and handles the translation from
 * c++ exceptions to Tcl errors.
 */
static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]);
  
void TclCommand::create_command(const string &command) {
  Tcl_CreateCommand(interp(), command.c_str(), reinterpret_cast<Tcl_CmdProc *>(TclCommand_wrapper), reinterpret_cast<ClientData>(this), NULL);    
}

static int TclCommand_wrapper (ClientData data, Tcl_Interp *interp, int argc, char *argv[]) {
  TclCommand *p = reinterpret_cast<TclCommand *>(data);

  if(argc > 1) {
    /* argv[0] is the command name */
    list<string> args(argv + 1, argv + argc);

    try {
        p->parse_from_string(args);
        
      if(!args.empty()) {
        std::stringstream ss;
        ss << "Unknown argument '" << args.front() << "' in '";

        for(int i = 0; i < argc; i++)
          ss << argv[i] << " ";
        ss << "'";
        throw std::runtime_error(ss.str());
      }
    } catch(std::exception &e) {
      Tcl_AppendResult(interp, e.what(), 0);
      return TCL_ERROR;
    }
  } else {
    Tcl_AppendResult(interp, p->print_to_string().c_str(), 0);
  }
  return TCL_OK;
}
