FROM ubuntu:focal

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y \
    gpg debmake debhelper devscripts equivs \
    distro-info-data distro-info

COPY entrypoint.sh /entrypoint.sh
COPY build.sh /build.sh

ENTRYPOINT ["/entrypoint.sh"]