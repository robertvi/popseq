#!/usr/bin/env python3

'''
create simple simulated one chromosome genomes for two diploid parents
and N simulated progeny
'''

import argparse
import random

class individual:
    def __init__(self,name,seq1,seq2):
        self.name = name
        self.chrm = [seq1,seq2]

    def __str__(self):
        return "Name: " + self.name + "\n" + self.chrm[0][:100] + "...\n" + self.chrm[1][:100] + '...'

    def get_chrm(self,i):
        return self.chrm[i]

def random_seq(size):
    'return a random sequence of the letters ATCG'
    return ''.join([random.choice("ATCG") for x in range(size)])

def generate_snp_list(genome_size,snps,debug=False):
    'return a list of positions and the two position variants at that location'

    #possible allele patterns for the two parents (0=varA,1=varB)
    options = ['0011','1100',
               '0100','0111','1011','1000',
               '0001','1101','1110','0010',
               '0101','1010','0110','1001']

    #determine positions of SNPs, ensure none coincide
    posn = [x for x in range(genome_size)]
    if not debug: random.shuffle(posn)
    posn = posn[:snps]

    snp_list = []
    for i in posn:
        #pick variantA at random (overwrites any exising base)
        varA = 'ATCG'[ random.choice([0,1,2,3]) ]

        #pick variantB being certain it doesn't match variantA
        indA = 'ATCG'.index(varA)
        adjB = random.choice([1,2,3])
        varB = 'ATCG'[(indA+adjB)%4]

        #decide which alleles each parent will have
        pattern = random.choice(options)

        #store position and the two variants
        snp_list.append([i,varA,varB,pattern])

    return snp_list

def generate_parents(seq,snp_list):
    'generate two parental sequences from a base sequence and snplist'

    p1 = [bytearray(seq,encoding='utf-8'),bytearray(seq,encoding='utf-8')]
    p2 = [bytearray(seq,encoding='utf-8'),bytearray(seq,encoding='utf-8')]

    #assign the alleles from each snp in the list
    for posn,varA,varB,pattern in snp_list:
        vars = varA+varB
        p1[0][posn] = ord(vars[int(pattern[0])])
        p1[1][posn] = ord(vars[int(pattern[1])])
        p2[0][posn] = ord(vars[int(pattern[2])])
        p2[1][posn] = ord(vars[int(pattern[3])])

    p1[0] = p1[0].decode('utf-8')
    p1[1] = p1[1].decode('utf-8')
    p2[0] = p2[0].decode('utf-8')
    p2[1] = p2[1].decode('utf-8')

    return individual("p1",p1[0],p1[1]),individual("p2",p2[0],p2[1])

def generate_test_parents(size):
    'generate parents of polyA,T,C,G for testing'

    return individual("p1",'A'*size,'T'*size),individual("p2",'C'*size,'G'*size)

def crossover(seq1,seq2,crossovers):
    'return a single sequence incorporating the requests number of crossovers'

    #pick positions of the crossovers
    size = len(seq1)
    posn = [x for x in range(1,size)] #[1,...size-1]
    random.shuffle(posn)
    posn = posn[:crossovers] #crossover before the offset position(s) given

    #begin as copy of both of the parent's chromosomes
    seqA = seq1[:]
    seqB = seq2[:]

    #perform each crossover
    for i in posn:
        temp = seqA[:]
        seqA = seqA[:i] + seqB[i:]
        seqB = seqB[:i] + temp[i:]

    #return one of the two chromosomes at random
    return random.choice([seqA,seqB])

def generate_progeny(name,p1,p2,crossovers):
    'generate genome sequence of progeny from the parents'

    seq1 = crossover(p1.get_chrm(0),p1.get_chrm(1),crossovers)
    seq2 = crossover(p2.get_chrm(0),p2.get_chrm(1),crossovers)

    return individual(name,seq1,seq2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('--progeny', type=int,required=True, help='Number of progeny')
    parser.add_argument('--snps', type=int,required=True, help='Number of SNPs')
    parser.add_argument('--genome-size', type=int, default=100000, help='Size of one chromosome in bases')
    parser.add_argument('--crossovers', type=int, default=2, help='Number of crossovers')

    args = parser.parse_args()

    print("progeny=%d"%args.progeny)
    print("genome_size=%d"%args.genome_size)
    print("snps=%d"%args.snps)

    #create random chromosome of required size
    base_seq = random_seq(args.genome_size)

    #create list of snp positions, their alleles and parental states
    snp_list = generate_snp_list(args.genome_size,args.snps)

    #generate parents' genome sequence
    parent1,parent2 = generate_parents(base_seq,snp_list)

    print(parent1)
    print(parent2)

    #generate progeny genome sequences
    progeny_list = []

    for i in range(args.progeny):
        progeny_list.append(generate_progeny("progeny%03d"%i,parent1,parent2,args.crossovers))
        print(progeny_list[-1])
