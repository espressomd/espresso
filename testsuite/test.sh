#!/bin/sh

####### guess srcdir from script location for manual execution
if test ! -f $srcdir/test.sh; then
    srcdir=`dirname $0`
fi

if test -x ../Espresso; then
    ESPRESSO=`pwd`/../Espresso
fi

# bail out if the environment is not correctly set
bail=
if test ! -f "$srcdir/test.sh"; then
    echo "please specify the variable srcdir pointing to the testsuite directory" 1>&2
    bail=y
fi
if test ! -x "$ESPRESSO"; then
    echo "please specify the variable ESPRESSO pointing to the Espresso program" 1>&2
    bail=y
fi
if test -z "$PROCESSORS"; then
    echo "please specify the processor numbers in PROCESSORS" 1>&2
    bail=y
fi
if test -z "$TESTCASES"; then
    echo "please specify the testcases to run in TESTCASES" 1>&2
    bail=y
fi
if test -n "$bail"; then
    exit -1
fi

# here go the error messages
errf=_test.sh_error.$$

# absolute path of source directory
srcpath=`cd $srcdir && pwd`
####### link all files from the src directory to the run directory, if they are not identical
if test `pwd` != $srcpath; then
    for f in $srcpath/*; do ln -f -s $f .; done
fi

# list of tests that are not supported by this version
blacklist=
# and what is missing for what test
missing=
for np in $PROCESSORS; do
    for f in $TESTCASES; do
	ignore=0
	for ft in $blacklist; do
	    if test "$ft" = "$f"; then
		ignore=1;
		break;
	    fi
	done
	# only not blacklisted tests
	if test $ignore -eq 0; then
	    # this is removed if the script runs through
	    echo "execution of script failed at unexpected point" > $errf
	    $ESPRESSO $srcpath/$f $np $errf
	    if test -f $errf; then
		if grep -q -e "^not compiled in:" $errf; then
		    missing="$missing$f: `cat $errf`\n"
		    blacklist="$blacklist $f"
		    rm -f $errf
		else
		    echo "Testcase $f failed for $np nodes with error `cat $errf`."
		    rm -f $errf
		    exit -666
		fi
	    fi
	else
	    echo "Test $f is blacklisted, ignoring..."
	fi
    done
done

if test "$missing" != ""; then
    echo -e "\n\n===============================================\n\n"
    echo -e "                Tests not done:\n                "
    echo -e $missing
fi

echo -e "\n\n===============================================\n\n"
echo -e "   Congratulations! ESPResSo seems to be ok."
echo -e "\n\n===============================================\n\n"
