if ((`ps -o nice= $$` < 5)); then 
    renice -n 5 $$
fi

function start() {
    echo "START $1"
}

function end() {
    echo "END $1"
}

function use_myconfig() {
    myconfig=testsuite/configs/myconfig-$1.h
    if [ ! -e "$myconfig" ]; then
        echo "$myconfig does not exist!"
        exit 1
    fi
    echo "Using $myconfig."
    cp $myconfig myconfig.h
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
    make dist dist-xz
    end "DIST"
}

function bootstrap() {
    start "BOOTSTRAP"
    ./bootstrap.sh
    end "BOOTSTRAP"
}
