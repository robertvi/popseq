#!/usr/bin/python

#
# group RGxHA lmxll markers into bins using the popseq progeny
#

import sys
#file matching popseq sample names to map progeny names in header of haplotype files
#sample_key = 'popseq2progenyname.csv' #all progeny
#sample_key = 'popseq2progenyname_best.csv' #excluding RGXH174a
sample_key = sys.argv[1]

#format string to match all haplotype files
conmap_fmt = '/home/vicker/octoploid_mapping/consensus_map4/popn_RGxHA/map/haplotypes/%s_haplotypes.csv'

#construct list of chromosome names
numb = '1234567'
letter = 'ABCD'
chrm_list=[]
for n in numb:
    for l in letter:
        chrm_list.append(n+l)

#load list of popseq progeny names
#must be in same order as the columns in the popseq kmer data
#labels must exactly match the columns in the map haplotype files
popseq_prog = {}
ctr = 0
f = open(sample_key)
for line in f:
    tok = line.strip().split(',')
    popseq_prog[tok[1]] = [ctr,0] #hash of which progeny are to be included and their position in the order
    ctr += 1
f.close()

code2info = {}
for chrm in chrm_list:
    #print chrm
    f = open(conmap_fmt%chrm)
    header = f.readline().strip().split(',')

    #record column offset of each popseq progeny maternal haplotype in the map file
    for i in xrange(len(header)):
        if header[i] in popseq_prog:
            popseq_prog[header[i]][1] = i

    for line in f:
        tok = line.strip().split(',')

        #process only lmxll markers
        if tok[2] != '<lmxll>': continue

        cm = float(tok[1])

        phase = tok[3][1]  # from {0-} or {1-}

        #create binmap barcode, 0=hom 1=het
        raw0 = [None] * len(popseq_prog)
        raw1 = [None] * len(popseq_prog)
        for x in popseq_prog:
            i,j = popseq_prog[x] #required position in binmapcode, column offset in map file
            mat_call = tok[j]
            pat_call = tok[j+1]
            assert mat_call in 'AB' and pat_call in 'AB'
            if mat_call == pat_call: #hom
                raw0[i] = '0' #raw calls encoded as 0=hom (lacking unique kmer) 1=het (containing unique kmer)
                raw1[i] = '1' #raw calls encoded as 1=hom 0=het
            else: #het
                raw0[i] = '1' #raw calls encoded as 0=hom 1=het
                raw1[i] = '0' #raw calls encoded as 1=hom 0=het

        assert not None in raw0

        #"canonical" binmap code encoded wrt progeny1
        #record phase as if progeny1 == 0
        if raw0[0] == '0':
            code = ''.join(raw0[1:])
        else:
            code = ''.join(raw1[1:])
            phase = '0' if phase == '1' else '1' #flip phase

        #initialise new record if required
        if not code in code2info: code2info[code] = {}
        if not chrm in code2info[code]: code2info[code][chrm] = [0,0,phase,cm,cm,0.0] #0-count,1-count,0-phase,min_cm,max_cm,total_cm

        #increment appropriate counter
        if raw0[0] == '0': code2info[code][chrm][0] += 1
        else:              code2info[code][chrm][1] += 1

        #record cm position
        if cm < code2info[code][chrm][3]: code2info[code][chrm][3] = cm
        if cm > code2info[code][chrm][4]: code2info[code][chrm][4] = cm
        code2info[code][chrm][5] += cm

    f.close()

if False:
    #display total number of of codes
    print 'codes',len(code2info)
    #find ambiguous codes
    for code in code2info:
        if len(code2info[code]) > 1:
            print code
            for x in code2info[code]:
                print x,code2info[code][x]
else:
    #output: code chrm count0 count1 phase0 mincm maxcm meancm
    for code in code2info:
        for x in code2info[code]:
            print code,x,

            #convert total cm into mean cm
            code2info[code][x][5] = code2info[code][x][5] / float(code2info[code][x][0]+code2info[code][x][1])
            for y in code2info[code][x]: print y,
            print

