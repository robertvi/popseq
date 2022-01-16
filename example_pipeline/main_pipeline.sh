#!/usr/bin/env bash

#
# example pipeline - a complete list of all the steps of the pipeline
# in the order they need to be executed
# note: you may not want to simply run this script
# but rather execute each step one at a time
# using any appropriate HPC Job Scheduler commands etc
#
# copy this script into a suitable folder and modify the key path variables below
# set them up to point to your raw input, working data and script folders
#

set -eu

# ==============================================================
# modify these paths to point to your:
#  - raw input data folder (DATADIR)
#  - working/scratch data folder (WORKDIR)
#  - popseq installation folder (INSTALLDIR)
#  - dsk installation folder (DSKDIR)
# ==============================================================

export DATADIR=/home/vicker/bioinformatics/raw_data
export WORKDIR=/home/vicker/bioinformatics/working_data
export INSTALLDIR=/home/vicker/bioinformatics/software/popseq
export DSKDIR=/home/vicker/bioinformatics/software/dsk

export PATH=${INSTALLDIR}/scripts:${PATH}
export PATH=${INSTALLDIR}/build/bin:${PATH}
export PATH=${DSKDIR}/build/bin:${PATH}

cd ${WORKDIR}
pwd
