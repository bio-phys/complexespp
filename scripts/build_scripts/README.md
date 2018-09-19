This folder contains the build scripts used to generate a Docker image for
continues integration. But these scripts can also be used for convenience to
install complexes and it's dependencies.

# Install Complexes
To have a normal build run the following commands

    ./build_deps.sh
    ./build_dynamic.sh

This will install complexes in `~/complexes-pp/bin`. To change this path modify
`build.cfg`. To build a static binary replace the last command with
`./build_static.sh`.

# Build Docker Image

To create the docker image run `sudo ./build-docker-image.sh`. This commands
needs to be run with sudo rights because docker doesn't work for normal users.

## Test Docker Image locally

You can test the Docker image by creating a short scripts that downloads the
complexes repo and runs the tests.

    git clone https://username:password@gitlab.mpcdf.mpg.de/MPIBP-Hummer/complexes-pp.git /builds/gitlab-org/complexes-pp
    cd /builds/gitlab-org/complexes-pp
    # copy of script lines in `.gitlab-ci.yml` from root folder
    bash scripts/gitlab-ci/gitlab-python-ci.sh
    bash scripts/gitlab-ci/gitlab-complexes-ci.sh
    bash scripts/gitlab-ci/gitlab-tutorials-ci.sh

In the above script replace `username` and `password`. This file will have your
password in clear text. So make sure that only you have reading rights for and
and NEVER EVER add it to the git repository.

To run the build as a test use.

    sudo docker run -i --name build_test complexes-builder /bin/bash < build_script
