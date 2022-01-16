#!/bin/bash

#
# popseq analysis of rgxha progeny
# 20 progeny were whole genome sequenced
#

#
# raw barcodes are the kmer presence/absence in each of the n progeny
# encoded as 1 or 0 respectively
# canonical barcodes are the last n-1 values encoded with respect to the
# value of the first progeny's value, 0=same as first progeny 1=different from
# this is so that codes will match regardless of phase
# (see popseq_utils.cc calls2code function)
# eg if there are n=5 progeny:
# raw code 01100 becomes 1100
# raw code 10011 becomes 1100
# raw code 11100 becomes 0011
# the raw code of the first progeny is retained as well to allow translation back


set -eu

export DATADIR=/data/seq_data/external/20180305_novogene_octoseq_popseq
export OCTODIR=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq
#export WORKDIR=/home/vicker/octoseq_project/popseq
export WORKDIR=/data/scratch/vicker/octoseq_project/popseq

export PATH=/home/vicker/git_repos/popseq:${PATH}
export PATH=/home/vicker/programs/dsk-2.1.0-Linux/bin:${PATH}

cd ${WORKDIR}

#=======================================================================
# count reads per sample
#=======================================================================
ls -1 ${DATADIR} > sample_list
rm -f linecounts

#count lines in fastq files
for x in $(cat sample_list)
do
    for y in $(ls -1 ${DATADIR}/${x}/*.fq.gz)
    do
        echo ${y}
        echo -n "${x} " >> linecounts
        zcat ${y} | wc --lines >> linecounts
    done
done

#convert lines into total bases
linect2bases.py > base_counts

#=======================================================================
# count kmers in RG and Hapil
#=======================================================================
#make list of Redgauntlet Illumina paired end fastq files
ls -1 ${OCTODIR}/PE/*.gz > rg_fastq_list

#make list of Hapil Illumina paired end fastq files
ls -1 ${OCTODIR}/cultivar_seq/Hapil*.gz > ha_fastq_list

for x in ha rg
do
    #run dsk
    dsk -file ${x}_fastq_list -kmer-size 31 -abundance-min 1 -max-memory 5000

    #extract histogram counts
    h5dump -y -d histogram/histogram ${x}_fastq_list.h5 > ${x}_histo_k31
    gzip ${x}_fastq_list.h5 #ha=88G uncompressed rg  466G uncompressed

    #convert hdf5 output into tsv format
    cat ${x}_histo_k31 | grep "^\ *[0-9]" | tr -d " ," | paste -d, - - > ${x}_histo_k31.csv
done

#=======================================================================
# count kmers in progeny
#=======================================================================
for x in $(cat sample_list)
do
    #make list of paired end fastq files
    ls -1 ${DATADIR}/${x}/*.fq.gz > ${x}_list

    #run dsk
    dsk -file ${x}_list -kmer-size 31 -abundance-min 1 -max-memory 5000

    #extract histogram counts
    h5dump -y -d histogram/histogram ${x}_list.h5 > ${x}_histo_k31
    gzip ${x}_fastq_list.h5 #ha=88G uncompressed

    #convert hdf5 output into tsv format
    cat ${x}_histo_k31 | grep "^\ *[0-9]" | tr -d " ," | paste -d, - - > ${x}_histo_k31.csv
done

#=======================================================================
# plot kmer frequency histograms, identify first peak
#=======================================================================
plot_kmer_specs.R

#=======================================================================
# dump from h5 to text format, filter out most error kmers
#=======================================================================
#progeny
filter_progeny_kmers.sh

#hapil
dsk2ascii_stdout -file ha_fastq_list.h5 -out dummy \
    2> /dev/null \
    | grep -v ':' \
    | awk '$2 > 4' \
    | lzop \
    > ha_kmers.lzo

#redgauntlet
dsk2ascii_stdout -file rg_fastq_list.h5 -out dummy \
    2> /dev/null \
    | grep -v ':' \
    | awk '$2 > 20' \
    | lzop \
    > rg_kmers.lzo

#=======================================================================
# call genotypes from kmer counts
#=======================================================================
qsub run_process_kmers.sh #outputs bestdata_calls.lzo

#=======================================================================
# aggregate kmers into binmap codes, taking account of repulsion phase
#=======================================================================
#process_calls alldata_calls.lzo > all_code_counts

#codes 29201
process_calls bestdata_calls.lzo > best_code_counts
plot_code_counts.R

#=======================================================================
# group i90 RGxHA lmxll markers into bins
#=======================================================================
#alldata
#364 bins with no collisions

#351 bincodes with no collisions
rgxha_lmxll_bincodes.py popseq2progenyname_best.csv > best_lmxll_code_counts
plot_lmxll_codecounts.R

#=======================================================================
# match up popseq kmer bins with map bins
#=======================================================================
#alldata - all 20 progeny
#175 of 364 bins have exact matches
#3496558 kmers matched out of 99546949 1-RG 0-HG kmers

#bestdata - 19 best progeny
#165 bins of 351 have exact matches
#6,593,629 kmers matched out of 92,980,970
#assuming 6.6M exact matches out of 93M = 0.07 error free codes
#implies 0.13 error rate per call (0.13^19 = 6.6/93)
#implies 7% 0 error, 20% 1 error, 27% 2 error
matchup_bincodes.py best_lmxll_code_counts best_code_counts

#matchup and output consolidated info into single file
#still using only exact matches
lzop -dcf bestdata_calls.lzo | final_matchup_bincodes.py best_lmxll_code_counts | lzop > best_consolidated.lzo

#count kmers per bin ready for R plotting
#avoids using geom_histogram to do the counting which gets confused by zero counts when also log scaling
count_final_matchups.py best_consolidated > best_matchup_counts

plot_kmer2bins.R

#=======================================================================
# fuzzy matchup of popseq kmer bins with map bins in one step
#=======================================================================
#56089334 matches
fuzzy_match best_lmxll_code_counts bestdata_calls.lzo 2> /dev/null | lzop > best_fuzzy.lzo

mkfifo fuzzy_fifo
lzop -dcf best_fuzzy.lzo > fuzzy_fifo &
count_final_matchups.py fuzzy_fifo > best_fuzzy_counts
