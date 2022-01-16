#!/usr/bin/env python3

'''
flatten the sequences from a fastq file
so that each sequence is on a single line
'''

import sys

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: flatten_fastq.py <input_filename>")
    print("results are printed to stdout")
    exit()

input_file = sys.argv[1]

with open(input_file) as f:
    for line_count,line_data in enumerate(f):
        #print identifier
        if line_count%4 == 0:
            print(line_data[1:].strip(),end=' ')
        #print sequence
        elif line_count%4 == 1:
            print(line_data[1:].strip(),end=' ')
        #ignore second identifier
        elif line_count%4 == 2:
            #print(line_data[1:].strip(),end='\n')
            continue
        #print quality
        elif line_count%4 == 3:
            print(line_data[1:].strip(),end='\n')
