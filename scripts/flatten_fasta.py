#!/usr/bin/env python3

'''
flatten the sequences from a fasta file
so that each sequence is on a single line
'''

import sys

#provide usage information
if sys.argv[1] in ['-h','-H','--help','-?']:
    print("Usage: flatten_fasta.py <input_filename>")
    print("results are printed to stdout")
    exit()

input_file = sys.argv[1]
current_identifier = None
current_sequence = ''

with open(input_file) as f:
    for line_data in f:
        #beginning of a new record
        if line_data.startswith('>'):
            #output previous record if any
            if current_identifier != None:
                print(current_identifier+' '+current_sequence)

            #start new record
            current_identifier = line_data[1:].strip()
            current_sequence = ''
            continue

        #append sequence to current record
        current_sequence += line_data.strip()

#print out the final record
print(current_identifier+' '+current_sequence)
