#!/bin/sh
version=`awk '/^\* / { for(i=2;i<NF-1;i++) printf("%s ",$i); print $(NF-1); exit}' RELEASE_NOTES`
# protect &'s in some code names chosen
version=`echo $version | sed 's/&/\\\&/g'`
date=`awk '/[(][a-z]*[)]/ { for(i=2;i<NF;i++) printf("%s ",$i); gsub("[.]", "", $NF); print $NF; exit}' RELEASE_NOTES`
echo "generating configure.ac for Espresso $version, dated $date"
sed -e "s/@VERSION@/$version/" -e "s/@TIME_STAMP@/$date/" configure.ac.raw > configure.ac
echo "generating configure"
aclocal --acdir=config
autoconf
