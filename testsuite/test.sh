#!/bin/sh

for f in *.tcl; do
    if test -x $f; then
	if ! $f; then
	    set error 1
	    break
	fi
    fi
done

if test "$error" == "1"; then
    echo "Testcase $f failed."
fi