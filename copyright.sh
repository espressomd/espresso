#!/bin/tcsh
foreach i ( `grep "Copyright (c) 2002-2003" Makefile* internal/*.tcl samples/*.tcl scripts/*.tcl testsuite/*.tcl | awk '{print $1}' | sed 's|:#||'` )
    cat $i | sed 's/Copyright (c) 2002-2003/Copyright (c) 2002-2004/' > $i.new
    touch --reference=$i $i.new
    mv -f $i.new $i
end
foreach i ( `grep "Copyright (c) 2002-2003" *.[ch] | awk '{print $1}' | sed 's|://||'` )
    cat $i | sed 's/Copyright (c) 2002-2003/Copyright (c) 2002-2004/' > $i.new
    touch --reference=$i $i.new
    mv -f $i.new $i
end



