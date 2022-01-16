#!/usr/bin/env bash

#
# test the top level script
#

set -eu

parental_reads=10000
progeny_reads=1000

../scripts/create_test_genomes.py --progeny 10 --genome-size 1000 --snps 5

../scripts/generate_reads.py --flata parent1.flata --reads ${parental_reads} > parent1.flatq
../scripts/generate_reads.py --flata parent2.flata --reads ${parental_reads} > parent2.flatq

ls -1 progeny*.flata > progeny_list
for x in $(cat ./progeny_list)
do
    ../scripts/generate_reads.py --flata ${x} --reads ${progeny_reads} > ${x/\.flata/}.flatq
done

for x in *.flatq
do
    ../scripts/unflatten_fastq.py ${x} > ${x/\.flatq/.fastq}
done
