#!/usr/bin/env bash

#
# example pipeline - a complete list of all the steps of the pipeline
# in the order they need to be executed
# note: you may not want to simply run this script
# but rather execute each step one at a time
# using the appropriate HPC Job Scheduler commands
#
# copy this script into a suitable folder and modify the key path variables below
# set them up to point to your raw input, working data and script folders
#

set -eu

# ==============================================================
# modify these paths to point to your input data, working data
# and installation folders
# ==============================================================

# export DATADIR=/data/seq_data/external/20180305_novogene_octoseq_popseq
# export OCTODIR=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq
# export WORKDIR=/data/scratch/vicker/octoseq_project/popseq
# export PATH=/home/vicker/git_repos/popseq:${PATH}
# export PATH=/home/vicker/programs/dsk-2.1.0-Linux/bin:${PATH}

export DATADIR=/data/seq_data/external/20180305_novogene_octoseq_popseq
export OCTODIR=/home/groups/harrisonlab/project_files/fragaria_x_ananassa/octoseq
export WORKDIR=/data/scratch/vicker/octoseq_project/popseq
export PATH=/home/vicker/git_repos/popseq:${PATH}
export PATH=/home/vicker/programs/dsk-2.1.0-Linux/bin:${PATH}

cd ${WORKDIR}
