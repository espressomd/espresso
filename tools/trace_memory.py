#!/usr/bin/python
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2008 Axel Arnold
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import print_function
import sys
import re

maxbacktrace=5

f=open(sys.argv[1], "r")
if len(sys.argv) > 2:
    n=int(sys.argv[2])
else:
    n=0

# regular expressions
re_start = re.compile(r"^%d: (?P<op>[a-z]+) (?P<args>.*)" % n)

allocated = {}

linenr=0
for line in f:
    linenr = linenr + 1
    if linenr % 1000 == 0:
        sys.stderr.write(".")
    match = re_start.match(line)
    if match == None: continue
    op = match.group('op')
    args = match.group('args').split(" ")
    if op == "alloc":
        size = args[0]
        addr = args[2]
        src  = [args[4]]
        allocated[addr] = (size, src)
    elif op == "realloc":
        old  = args[0]
        addr = args[2]
        size = args[4]
        src  = [args[6]]
        if old == "(nil)":
            pass
        elif old in addr:
            prev = allocated[old][1][:maxbacktrace-1]
            src.extend(prev)
            del allocated[old]
        else:
            src.append("unmanaged source " + old)
        allocated[addr] = (size, src)
    elif op == "free":
        addr = args[0]
        src  = args[2]
        if addr == "(nil)":
            pass
        elif addr in allocated:
            del allocated[addr]
        else:
            print(("\n" + addr + " freed at " + src + ", but never allocated\n"))

print("\n")
for (addr,info) in list(allocated.items()):
    s = info[0] + " @ " +  addr + " allocated at " + info[1][0]
    for loc in info[1][1:]:
        s += ", from " + loc
    print(s)
