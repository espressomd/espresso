#!/bin/sh


errf=_test.sh_error.$$

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
	    cdir=`pwd`
	    odir=$cdir/..
	    ( cd $srcdir; $odir/Espresso $f $np $cdir/$errf )
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
