#!/bin/sh

for np in 1 2 3 4 6 8; do
    for f in *.tcl; do
	if test -x $f; then
	    if ! $f $np; then
		set error 1
		break
	    fi
	fi
    done
done

if test "$error" == "1"; then
    echo "Testcase $f failed."
fi