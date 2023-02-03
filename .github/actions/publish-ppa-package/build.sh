#!/bin/bash

set -o errexit -o pipefail -o nounset

REPOSITORY=$INPUT_REPOSITORY
GPG_PRIVATE_KEY="$INPUT_GPG_PRIVATE_KEY"
GPG_PASSPHRASE=$INPUT_GPG_PASSPHRASE
PKGDIR=$INPUT_PKGDIR
IS_NATIVE=$INPUT_IS_NATIVE
SERIES=$INPUT_SERIES

assert_non_empty() {
    name=$1
    value=$2
    if [[ -z "$value" ]]; then
        echo "::error::Invalid Value: $name is empty." >&2
        exit 1
    fi
}

assert_non_empty inputs.repository "$REPOSITORY"
assert_non_empty inputs.gpg_private_key "$GPG_PRIVATE_KEY"
assert_non_empty inputs.gpg_passphrase "$GPG_PASSPHRASE"
assert_non_empty inputs.pkgdir "$PKGDIR"

if [[ -z "$SERIES" ]]; then
    SERIES=$(distro-info --supported)
fi

echo "::group::Importing GPG private key..."
GPG_KEY_ID=$(echo "$GPG_PRIVATE_KEY" | gpg --import-options show-only --import | sed -n '2s/^\s*//p')
echo $GPG_KEY_ID
echo "$GPG_PRIVATE_KEY" | gpg --batch --passphrase "$GPG_PASSPHRASE" --import
echo "::endgroup::"

for s in $SERIES; do
    ubuntu_version=$(distro-info --series $s -r)
    version="${ubuntu_version:0:5}"

    echo "::group::Building deb for: $version ($s)"

    cd $PKGDIR
    if [[ -z "$IS_NATIVE" ]]; then
        format_file="debian/source/format"
        if [[ ! -f "$format_file" ]] || grep -q "native" "$format_file"; then
            IS_NATIVE="true"
        fi
    fi

    if [[ "$IS_NATIVE" == "true" ]]; then
        echo "Making native package..."
        debmake -n
    else
        echo "Making non-native package..."
        debmake
    fi

    mk-build-deps --install --remove --tool='apt-get -o Debug::pkgProblemResolver=yes --no-install-recommends --yes' debian/control

    sed -i"" -re "s/\s\w+;/ $s;/" -re "s/\((.+)-.+\)/(\1-ppa1~ubuntu${version})/" debian/changelog 

    if [[ "$IS_NATIVE" == "true" ]]; then
        echo "y" | debuild -S -sa \
            -k"$GPG_KEY_ID" \
            -p"gpg --batch --passphrase "$GPG_PASSPHRASE" --pinentry-mode loopback"
    else
        debuild -S -sa \
            -k"$GPG_KEY_ID" \
            -p"gpg --batch --passphrase "$GPG_PASSPHRASE" --pinentry-mode loopback"
    fi
    dput $REPOSITORY ../*.changes

    rm -rf ../*.{changes,build,buildinfo,deb,ddeb,dsc}
    echo "::endgroup::"
done