#!/usr/bin/env python3

'''
create simple simulated one chromosome genomes for two diploid parents
and N simulated progeny
'''

import argparse
import random

def random_seq(size):
    'return a random sequence of the letters ATCG'
    return ''.join([random.choice("ATCG") for x in range(size)])

def apply_snps(seq,nsnps):
    'return the input string with nsnps mutated to a different base'

    #determine positions to apply SNPs
    posn = [x for x in range(len(seq))]
    #random.shuffle(posn)
    posn = posn[:nsnps]

    mutant = bytearray(seq,encoding='utf-8')

    for i in posn:
        #find current base
        wt = 'ATCG'.index(chr(mutant[i]))

        #adjust current base by 1,2 or 3 places in the list ATCG with wrap
        adj = random.choice([1,2,3])

        mutant[i] = ord('ATCG'[(wt+adj)%4])
    return mutant.decode('utf-8')

class individual:
    def __init__(self,name,seq1,seq2):
        self.name = name
        self.chrm1 = seq1[:]
        self.chrm2 = seq2[:]

    def __str__(self):
        return "Name: " + self.name + "\n" + self.chrm1[:80] + "...\n" + self.chrm2[:80] + '...'

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--progeny', type=int,required=True, help='Number of progeny')
    parser.add_argument('--snps', type=int,required=True, help='Number of SNPs')
    parser.add_argument('--genome-size', type=int, default=100000, help='Size of one chromosome in bases')

    args = parser.parse_args()

    print("progeny=%d"%args.progeny)
    print("genome_size=%d"%args.genome_size)
    print("snps=%d"%args.snps)

    #create random chromosome of required size
    base_seq = random_seq(args.genome_size)
    snps_seq = apply_snps(base_seq,args.snps)

    parent1 = individual("p1",base_seq,snps_seq)
    print(parent1)
