#!/bin/sh
TESTCASES="lj.tcl"
# lj-cos.tcl FENE.tcl harmonic.tcl... constraints thermosim energy pressure rotation gay-berne

for np in 1 2 3 4 6 8; do
    for f in $TESTCASES; do
	if ! $f $np; then
	    set error 1
	    break
	fi
    done
done

if test "$error" == "1"; then
    echo "Testcase $f failed for $np nodes."
fi