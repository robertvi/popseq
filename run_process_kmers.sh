#$ -S /bin/bash
#$ -l h_vmem=20G
#$ -l mem_free=20G
#$ -l virtual_free=20G
#$ -l h_rt=9999:00:00
#$ -l h=blacklace08.blacklace
#$ -cwd
###$ -pe smp 3
###$ -t 1-20

#
# convert from kmer counts into genotype calls
# used about 10GB RAM

export PATH=/home/vicker/git_repos/popseq:${PATH}

export WORKDIR=/data/scratch/vicker/octoseq_project/popseq

cd ${WORKDIR}

#original calling process for progeny using progeny_peak_fraction parameter
#process_kmers subsample0.01_config | lzop > subsample0.01_calls.lzo
#process_kmers all_config | lzop > alldata_calls.lzo

#simplified calling process ignoring progeny_peak_fraction parameter
#where 3 indicates missing calls which should be excluded
#also excluding the worst progeny RGXH174a
process_kmers best_config | lzop > bestdata_calls.lzo
