#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y atlas-devel
yum install -y zlib-devel

# Compile wheels
for PYBIN in /opt/python/cp3[56]*/bin; do
    "${PYBIN}/pip" install -r /io/dev-requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/poreplex*.whl; do
    auditwheel repair "$whl" -w /io/dists/
done
