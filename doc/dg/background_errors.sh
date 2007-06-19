#!/bin/sh
# Parse the background-error lines from the source files

AWK_SCRIPT=$SRCDIR/background_errors.awk
UNSORTED=background_errors.unsorted
SORTED=background_errors.sorted
DOC=background_errors.doc

$AWK -f $AWK_SCRIPT "$@" > $UNSORTED
sort $UNSORTED > $SORTED

# OUTPUT
cat <<EOF > $DOC
/** \\page background_errors background_errors resolved
<ul>
EOF

cat $SORTED >> $DOC

cat <<EOF >> $DOC
</ul>
*/
EOF
