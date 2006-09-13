#!/bin/sh -x
aclocal -I config
autoheader
automake --add-missing --force-missing --copy
autoconf
