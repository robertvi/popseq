#$ -S /bin/bash
#$ -l h_vmem=12G
#$ -l mem_free=12G
#$ -l virtual_free=12G
#$ -l h_rt=9999:00:00
###$ -l h=blacklace02.blacklace
#$ -cwd
###$ -pe smp 3
###$ -t 1-20

export DATADIR=/data/seq_data/external/20180305_novogene_octoseq_popseq
export OCTODIR=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq
#export WORKDIR=/home/vicker/octoseq_project/popseq
export WORKDIR=/data/scratch/vicker/octoseq_project/popseq/rg_recount

export PATH=/home/vicker/git_repos/popseq:${PATH}
export PATH=/home/vicker/programs/dsk-2.1.0-Linux/bin:${PATH}
#export PATH=/home/vicker/programs/fastx_toolkit-0.0.14/src/fastx_quality_stats:${PATH}

cd ${WORKDIR}

#run dsk
dsk -file rg_fastq_list -kmer-size 31 -abundance-min 1 -max-memory 5000

#extract histogram counts
h5dump -y -d histogram/histogram rg_fastq_list.h5 > rg_histo_k31

#convert hdf5 output into tsv format
cat rg_histo_k31 | grep "^\ *[0-9]" | tr -d " ," | paste -d, - - > rg_histo_k31.csv

#compress
gzip ${x}_fastq_list.h5

