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

def generate_snp_list(genome_size,snps):
    'return a list of positions and the two position variants at that location'

    #determine positions of SNPs, ensure none coincide
    posn = [x for x in range(genome_size)]
    random.shuffle(posn)
    posn = posn[:snps]

    snp_list = []
    for i in posn:
        #pick variantA at random (overwrites any exising base)
        varA = 'ATCG'[ random.choice([0,1,2,3]) ]

        #pick variantB being certain it doesn't match variantA
        indA = 'ATCG'.index(varA)
        adjB = random.choice([1,2,3])
        varB = 'ATCG'[(indA+adjB)%4]

        #store position and the two variants
        snp_list.append([i,varA,varB])

def generate_parents(seq,snp_list):
    'generate two parental sequences from a base sequence and snplist'
    return individual("p1",seq,seq),individual("p2",seq,seq)

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
    snp_list = generate_snp_list(args.genome_size,args.snps)

    parent1,parent2 = generate_parents(base_seq,snp_list)
    print(parent1)
    print(parent2)
