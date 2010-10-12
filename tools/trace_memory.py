#!/usr/bin/python
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
            print "\n" + addr + " freed at " + src + ", but never allocated\n"

print "\n"
for (addr,info) in allocated.iteritems():
    s = info[0] + " @ " +  addr + " allocated at " + info[1][0]
    for loc in info[1][1:]:
        s += ", from " + loc
    print s
