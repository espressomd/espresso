FROM ubuntu:focal

ENV DEBIAN_FRONTEND noninteractive
run apt update
run apt install -y dpkg-dev build-essential

COPY entrypoint.sh /entrypoint.sh
