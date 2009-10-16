#!/bin/tcsh

# foreach i ( `grep "Copyright (c) 2002-2003" Makefile* internal/*.tcl samples/*.tcl scripts/*.tcl testsuite/*.tcl | awk '{print $1}' | sed 's|:#||'` )
#     cat $i | sed 's/Copyright (c) 2002-2003/Copyright (c) 2002-2004/' > $i.new
#     touch --reference=$i $i.new
#     mv -f $i.new $i
# end

# update the copyright of all .c and .h files
foreach i ( `grep "Copyright (c) 2002-2006" *.[ch] | awk '{print $1}' | sed 's|://||'` )
    echo $i
    cat $i | sed 's/Copyright (c) 2002-2006/Copyright (c) 2002-2009/' > $i.new
    touch --reference=$i $i.new
    mv -f $i.new $i
end



