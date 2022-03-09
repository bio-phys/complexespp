#!/bin/bash

# needs to be run as root for docker privileges, or use podman
docker build --tag=malinke/complexesbuilder:debian --file=Dockerfile.debian11 .
#docker build --tag=biophys/complexespp:builder-centos --file=Dockerfile.centos7 .
