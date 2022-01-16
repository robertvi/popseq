#!/usr/bin/python

#
# matchup the bincodes between RGxHA map and popseq kmers
#

import sys

mapfile = sys.argv[1]   #canonicalcode chrm 0-count 1-count 0-phase min_cm max_cm mean_cm
kmerfile = sys.argv[2]  #canonicalcode 0-count 1-count

code2info = {}

#load mapfile info
f = open(mapfile)
for line in f:
    tok = line.strip().split()
    
    code = tok[0]
    chrm = tok[1]
    count0 = int(tok[2])
    count1 = int(tok[3])
    phase0 = int(tok[4])
    mincm = float(tok[5])
    maxcm = float(tok[6])
    meancm = float(tok[7])
    
    code2info[code] = [chrm,count0,count1,phase0,mincm,maxcm,meancm]
    
f.close()

#print len(code2info)

matched = 0
total = 0
totaltotal = 0

f = open(kmerfile)
for line in f:
    tok = line.strip().split()
    
    code = tok[0]
    count0 = int(tok[1])
    count1 = int(tok[2])
    totaltotal += count0 + count1
    
    if code in code2info:
        matched += 1
        total += count0 + count1
f.close()

print matched
print totaltotal
print total
