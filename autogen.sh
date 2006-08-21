#!/bin/sh -x
aclocal -I config
autoheader
automake --foreign --add-missing
autoconf
