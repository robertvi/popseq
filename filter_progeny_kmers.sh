#!/bin/bash

export PATH=/home/vicker/git_repos/popseq:${PATH}
export PATH=/home/vicker/programs/dsk-2.1.0-Linux/bin:${PATH}

for x in $(cat sample_list)
do
    dsk2ascii_stdout -file ${x}_list.h5 -out dummy \
        2> /dev/null \
        | grep -v ':' \
        | awk '$2 > 2' \
        | lzop \
        > ${x}_kmers.lzo
done
