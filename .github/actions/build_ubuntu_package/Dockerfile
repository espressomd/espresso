FROM ubuntu:jammy

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y dpkg-dev build-essential

COPY entrypoint.sh /entrypoint.sh
