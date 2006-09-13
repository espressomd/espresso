#!/bin/sh -x
aclocal -I config
autoheader
automake -v --add-missing --copy
autoconf
