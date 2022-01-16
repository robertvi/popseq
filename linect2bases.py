#!/usr/bin/python

#
# convert line count data into total base pairs
#

linect = {}

with open("linecounts") as f:
    for line in f:
        sample,ct = line.strip().split()
        ct = int(ct)

        if not sample in linect: linect[sample] = 0

        linect[sample] += ct

print "samp bases coverage"
for sample in linect:
    linect[sample] /= 4    #convert from lines in fastq to (single) reads
    linect[sample] *= 150  #convert from reads to bases
    print sample,linect[sample],linect[sample]/800e6
