
#ifdef CPP_CONVERSIONS_H
#define CPP_CONVERSIONS_H

#include "tcl.h"
#include "Eigen/Core"

/* *
 * In this file we define a few generic functions to make returning 
 * complex c++ data easy. The idea is as follows: 
 *
 * We provide a set of functions 
 * void Tcl_Append(Tcl_Interp* interp, const T& what, bool is_in_list=0) 
 * that append what to the current result stored in the TCL interpreter.
 *
 * We provide a generic Tcl_Append function that uses the << operator
 * in combination with a stringstream. Any class that has an implentation
 * of <<(ostream) can be converted into a TCL result.
 * 
 * For certain data types we have special implementation: E.G. For 
 * double we use Tcl_PrintDouble, and a Vector3d from Eigen creates
 * a list with 3 entries. std::vector<T> is converted into a TCL list etc.
 *
 * The usage of lists is particularly interesting as this generic 
 * programming technique allows us also to use lists of lists, without
 * actually changing anything. One difficulty requires using a is_in_list
 * flag: TCL lists require having no braces {} when on "top level", but
 * lists in lists, need braces {}. Thus a simple list of list of integers 
 * looks like:
 * "{ 1 2 } { 3 4 }", without outer braces.
 *
 * Thus we pass the flag is_in_list to every Tcl_Append call, so that
 * inside the function is is possible to decide weather braces are necessary:
 * If is_in_list is set, outer braces are added, if the type itself
 * is a list.
 */

template <typename T>
void Tcl_Append(Tcl_Interp* interp, const T& what, bool is_in_list=0) {
  std::stringstream ss;
  ss << what;
  Tcl_AppendResult(interp, ss.str().c_str(), (char *)NULL);
}

void Tcl_Append(Tcl_Interp* interp, const double& what, bool is_in_list=0) {
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, what, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
}

template <size_t N>
void Tcl_Append(Tcl_Interp* interp, const char (&what)[N], bool is_in_list=0) {
  Tcl_AppendResult(interp, what, (char *)NULL);
}

void Tcl_Append(Tcl_Interp* interp, const Eigen::Vector3d &what, bool is_in_list=0) {
  if (is_in_list) 
    Tcl_Append(interp, "{");
  char buffer[TCL_DOUBLE_SPACE];
  Tcl_PrintDouble(interp, what[0], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, what[1], buffer);
  Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
  Tcl_PrintDouble(interp, what[2], buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  if (is_in_list) 
    Tcl_Append(interp, "}");
}

template <typename T>
void Tcl_Append(Tcl_Interp* interp, const std::vector<T> what, bool is_in_list=0) {
  append_iterable_to_tcl_result(interp, what, is_in_list);
}

template<class C>
void append_iterable_to_tcl_result(Tcl_Interp* interp, const C& what, bool is_in_list=0) {

  typename C::const_iterator it=what.begin();
  
  char* list = 0;
  char* position;
  
  if (is_in_list) 
    Tcl_Append(interp, "{");
  if (it!=what.end()) {
    Tcl_Append(interp, *it, true);
  }
  while (it!=what.end()) {
    Tcl_Append(interp, " ");
    Tcl_Append(interp, *it, true);
    ++it;
  }
  if (is_in_list) 
    Tcl_Append(interp, "}");
}


#endif
