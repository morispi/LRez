#!/bin/bash

set -e

make -j"$( nproc )" PREFIX="$( pwd )"
