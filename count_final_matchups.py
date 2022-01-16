#!/usr/bin/python

#
# count final matchups of kmers to bins
# output data ready to be plotted in R
# to avoid R getting confused with log10 scales of counts
#

import sys

#code chrm phase mincm maxcm meancm kmerseq rg ha p1 ... p19
inpfile = sys.argv[1]

count = {}               #count[chrm][phase][meancm] = count

f = open(inpfile)
for line in f:
    tok = line.strip().split()

    if tok[0] == 'code': continue #skip header

    if not tok[1] in count: count[tok[1]] = {} #chrm
    if not tok[2] in count[tok[1]]: count[tok[1]][tok[2]] = {} #phase

    if not tok[3] in count[tok[1]][tok[2]]:
        count[tok[1]][tok[2]][tok[5]] = 1
    else:
        count[tok[1]][tok[2]][tok[5]] += 1

f.close()

chrm_list = sorted(count.iterkeys())
for chrm in chrm_list:

    phase_list = sorted(count[chrm].iterkeys())
    for phase in phase_list:

        cm_list = sorted(count[chrm][phase].iterkeys(),key=float)
        for cm in cm_list:

            print chrm,phase,cm,count[chrm][phase][cm]
