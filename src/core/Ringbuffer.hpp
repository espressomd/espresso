/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#ifndef __RINGBUFFER_HPP
#define __RINGBUFFER_HPP

#include <deque>

template<typename T>
class Ringbuffer {
public:
  Ringbuffer(int _n) {
    if(_n <= 0)
      n = 1;
    else if (_n > d.max_size())
      n = d.max_size();
    else
      n = _n;
  };
  void push(T value) {
    if(d.size() >= n)
      d.pop_front();
    d.push_back(value);
  };
  typename std::deque<T>::iterator begin() { return d.begin(); };
  typename std::deque<T>::iterator end() { return d.end(); };
private:
  int n;
  std::deque<T> d;  
};

#endif
