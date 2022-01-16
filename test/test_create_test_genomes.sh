#!/usr/bin/env bash

#
# test the top level script
#

set -eu

../scripts/create_test_genomes.py --progeny 10 --genome-size 100 --snps 5
