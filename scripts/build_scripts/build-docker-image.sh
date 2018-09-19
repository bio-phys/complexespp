#!/bin/bash

# needs to be run as root for docker privileges
docker build --tag=biophys/complexespp:builder-debian --file=Dockerfile.debian9 .
docker build --tag=biophys/complexespp:builder-centos --file=Dockerfile.centos7 .
