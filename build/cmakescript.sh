#!/bin/bash

./cmakeclean.sh

cmake                                  \
  -DYAKL_CXX_FLAGS="${YAKL_CXX_FLAGS}" \
  -DYAKL_ARCH="${YAKL_ARCH}"           \
  -DLINK_FLAGS="${LINK_FLAGS}"         \
  ..


