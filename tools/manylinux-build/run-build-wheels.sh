#!/bin/sh
DOCKER_IMAGE=quay.io/pypa/manylinux1_x86_64
TOP_DIR=$(pwd)/../..

docker run --rm -v $TOP_DIR:/io $DOCKER_IMAGE /bin/bash /io/tools/manylinux-build/build-wheels.sh
