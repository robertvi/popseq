#!/usr/bin/env bash

#
# test the top level script
#

set -eu

../scripts/create_test_genomes.py --progeny x10 --genome-size 10000
