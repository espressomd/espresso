# Copyright (C) 2013 Olaf Lenz
#
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
if ((`ps -o nice= $$` < 5)); then 
    renice -n 5 $$
fi

function start() {
    echo "START $1"
}

function end() {
    echo "END $1"
}

function bootstrap() {
    cd $SRCDIR
    start "BOOTSTRAP"
    ./bootstrap.sh
    end "BOOTSTRAP"
    cd $BUILDDIR
}

function configure() {
    start "CONFIGURE"
    $SRCDIR/configure $@
    end "CONFIGURE"
}

function use_myconfig() {
    myconfig=$SRCDIR/testsuite/configs/myconfig-$1.hpp
    if [ ! -e "$myconfig" ]; then
        echo "$myconfig does not exist!"
        exit 1
    fi
    echo "Using $myconfig."
    cp $myconfig myconfig.hpp
}

function check() {
    start "TEST"
    # something should be done after ||, otherwise Jenkins will mark
    # job as failed
    make check || CHECK_UNSTABLE=1
    end "TEST"
}

function doc() {
    start "DOC"
    make doc
    end "DOC"
}

function dist() {
    start "DIST"
    make dist 
    make dist-xz
    end "DIST"
}

