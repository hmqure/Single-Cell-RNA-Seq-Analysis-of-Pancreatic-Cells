#!/bin/bash -l

#$ -P bf528
#$ -cwd
#$ -j y
#$ -pe mpi_16_tasks_per_node 16

# Setting up variable/environment
FASTQ_PATH=/projectnb/bf528/project_4_scrnaseq/fastq/
index=human_index
out=quant_o_test
sample=SRR3879604

module load salmon/1.1.0

#Use the preprocessed files as read 1 eg. SRRXXXXXXX_1_bc.fastq.gz

salmon quant -i $index -1 "${FASTQ_PATH}${sample}/${sample}_1_bc.fastq.gz" -2 "${FASTQ_PATH}${sample}/${sample}_2.fastq.gz"  --validateMappings -o $out

