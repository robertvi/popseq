#!/usr/bin/env python3

'''
convert flat sequences into fastq format
if no quality column insert a uniform placeholder quality
'''

import sys

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: unflatten_fastq.py <input_filename>")
    print("results are printed to stdout")
    exit()

input_file = sys.argv[1]

with open(input_file) as f:
    for line in f:
        column = line.strip().split()
        identifier = column[0]
        sequence = column[1]

        if len(column) >= 3:
            quality = column[2]
            assert(len(quality) == len(sequence))
        else:
            quality = '!' * len(sequence)

        print('@' + identifier)
        print(sequence)
        print('+' + identifier)
        print(quality)
