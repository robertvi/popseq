#!/usr/bin/env python3

'''
create simple simulated one chromosome genomes for two diploid parents
and N simulated progeny
'''

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--progeny', type=int,required=True, help='Number of progeny')
parser.add_argument('--genome-size', type=int, default=100000, help='Size of one chromosome in bases')

args = parser.parse_args()

print("progeny=%d"%args.progeny)
print("genome_size=%d"%args.genome_size)
