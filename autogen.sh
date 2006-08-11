#!/bin/sh
aclocal -I config
autoheader
automake --foreign --add-missing
autoconf
