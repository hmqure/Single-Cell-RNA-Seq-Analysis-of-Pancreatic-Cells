#!/bin/bash

#$ -P bf528
#$ -cwd
#$ -j y
#$ -pe mpi_16_tasks_per_node 16

module load salmon

transcripts=hg38_primary.fa

salmon index -t $transcripts -i human_index -k 31