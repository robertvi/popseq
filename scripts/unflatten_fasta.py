#!/usr/bin/env python3

'''
unflatten sequences into a fasta file
'''

import sys

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: unflatten_fasta.py <input_filename>")
    print("results are printed to stdout")
    exit()

input_file = sys.argv[1]
max_per_line = 80

with open(input_file) as f:
    for line in f:
        column = line.strip().split()
        identifier = column[0]
        sequence = column[1]

        print('>' + identifier)

        remaining = len(sequence)

        while remaining > 0:
            start = len(sequence) - remaining
            chunk_size = min(max_per_line,remaining)
            remaining -= chunk_size
            print(sequence[start:start+chunk_size])
