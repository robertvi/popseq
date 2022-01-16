#!/usr/bin/env python3

'''
generate simulated reads from sequences in flat format
'''

import sys
import random
from kmer_module import *
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--flata', required=True, help='Genome in flat fasta format')
    parser.add_argument('--reads', type=int,required=True, help='Number of reads')
    parser.add_argument('--length', type=int, default=100, help='Size of each read')
    parser.add_argument('--prob-error', type=float, default=0.01, help='Probability of read error per base')
    parser.add_argument('--prob-missing', type=float, default=0.01, help='Probability of uncallable base')
    parser.add_argument('--prob-revcomp', type=float, default=0.5, help='Probability of read coming from reverse compliment strand')

    args = parser.parse_args()


input_file = args.flata
reads = args.reads
length = args.length
prob_error = args.prob_error
prob_N = args.prob_missing
prob_revcomp = args.prob_revcomp

length_check_passed = False

#load in the source sequence(s)
sequence_list = []
with open(input_file) as f:
    for line in f:
        #split line into its columns
        column = line.strip().split()

        #for clarity store each column in a separate variable
        identifier = column[0]
        sequence = column[1]

        #we need at least one sequence longer than the read length
        if len(sequence) >= length:
            length_check_passed = True

        #append to list of sequences
        sequence_list.append([identifier,sequence])

assert(length_check_passed)

#generate the reads
for read_counter in range(reads):
    #pick a random sequence
    sequence_number = random.choice(range(len(sequence_list)))

    #pick a random starting position
    max_position = len(sequence_list[sequence_number][1]) - length
    start = random.randint(0,max_position)

    #extract the read
    read = sequence_list[sequence_number][1][start:start+length]

    #switch 50% to the reverse strand
    if random.random() < prob_revcomp:
        read = reverse_compliment(read)

    #split read into a list of bases for ease of modification
    base_list = [x for x in read]

    #apply errors and Ns
    for i,x in enumerate(base_list):
        if x == 'N': continue

        #convert to N
        if random.random() < prob_N:
            base_list[i] = 'N'

        #introduce a simple read error
        elif random.random() < prob_error:
            j = ('ATCG'.index(x) + random.randint(1,3)) % 4
            base_list[i] = 'ATCG'[j]

    ##convert back to a single string again
    read = ''.join(base_list)

    #print in flat format
    print(sequence_list[sequence_number][0]+'_'+str(read_counter)+' '+read)
