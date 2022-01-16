#!/usr/bin/python

#
# matchup the bincodes between RGxHA map and popseq kmers
# output kmer and all associated info for kmers which match a map bincode
#

import sys

mapfile = sys.argv[1]   #canonicalcode chrm 0-count 1-count 0-phase min_cm max_cm mean_cm
#kmerfile = sys.argv[2]  #kmerseq rg ha p1 ... p19, read from stdin

code2info = {}

#load mapfile info
f = open(mapfile)
for line in f:
    tok = line.strip().split()

    code = tok[0]
    chrm = tok[1]
    count0 = int(tok[2])
    count1 = int(tok[3])
    phase0 = int(tok[4]) #phase if progeny1==0
    mincm = float(tok[5])
    maxcm = float(tok[6])
    meancm = float(tok[7])

    #code2info[code] = [chrm,count0,count1,phase0,mincm,maxcm,meancm]
    code2info[code] = tok[1:]

f.close()

nprogeny = len(code) + 1

ct = 0
matched = 0
#f = open(kmerfile)
f = sys.stdin
header = f.readline().strip()

#header
print 'code chrm phase mincm maxcm meancm',header

for line in f:
    ct += 1
    #if ct %300000 == 0: print ct,matched

    tok = line.strip().split() #kmerseq rg ha p1 ... p19

    if tok[2] == '1': continue #ignore kmers present in hapil
    if '3' in tok: continue    #ignore kmers with missing calls

    #construct canonical code
    code = []
    for i in xrange(4,nprogeny+3): #n=19 progeny, column offsets 3+0,3+1...3+18
        if tok[i] == tok[3]: code.append('0')
        else:                code.append('1')

    code = ''.join(code)
    if not code in code2info: continue
    matched += 1

    #adjust phase to match status of progeny1
    info = code2info[code]

    chrm   = info[0]
    count0 = info[1]
    count1 = info[2]
    phase  = info[3]
    mincm  = info[4]
    maxcm  = info[5]
    meancm = info[6]

    progeny1 = tok[3]
    if progeny1 == '1': phase = '1' if phase=='0' else '0' #flip phase if progeny1 == '1'

    print code,               #map bincode
    print chrm,               #chromosome
    print phase,              #phase of kmer
    print mincm,maxcm,meancm, #cm position of bin
    print ' '.join(tok)       #kmerseq rg ha p1... p19

f.close()

#print ct,matched
