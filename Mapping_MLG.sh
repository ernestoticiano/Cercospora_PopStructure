#!/bin/bash
#PBS -l nodes=1:ppn=16
#PBS -l walltime=12:00:00
#PBS -N cflagellarisgenomes
#PBS -o mapping_MLG_stdout
#PBS -e mapping_MLG_stderr
#PBS -q single
#PBS -A hpc_richards04
#PBS -m e
#PBS -M edasilva@agcenter.lsu.edu

module load bwa
module load samtools
module load gatk
module load jdk
module load gnuparallel

cd /work/ets/GenomeResequencing

### Index reference genome using BWA; this can sometimes throw errors depending on where things are placed. ###
### I suggest moving the reference genome to the working directory (in the directory of the path above) ### 

bwa index gene.ref.fa

###Setup for loop using forward read FASTQ file###

for seq in /work/ets/GenomeResequencing/Trimmed_Reads/*_r1_trimmed.fq.gz
do

###Use sed to remove suffix and keep sample name###

name=$(echo $seq | sed 's/_r1_trimmed.fq.gz//')
echo $name

###Map reads to reference genome, adding read group information###

bwa mem -t 16 -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA\tLB:run1" gene.ref.fa ${name}_r1_trimmed.fq.gz ${name}_r2_trimmed.fq.gz > ${name}.sam

###Convert SAM to BAM format, sort, mark duplicate reads, and idex###
samtools view -@ 12 -b -S ${name}.sam > ${name}.bam
samtools sort -@ 12 ${name}.bam > ${name}.sorted.bam
gatk MarkDuplicates -I ${name}.sorted.bam -O ${name}.sorted.md.bam -M ${name}_marked_dup_metrics.txt
samtools index ${name}.sorted.md.bam
done

### Make sure you're left with {name}.sorted.md.bam files. ###
