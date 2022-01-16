#!/usr/bin/env bash

#
# run all tests, exit with non-zero code on the first failure
# run from within the test folder

set -eu

trap 'echo "Errorcode $? on line $LINENO" ; exit 1' ERR

#python script tests
./test_create_test_genomes.sh
