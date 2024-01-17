#!/bin/sh

#SBATCH --job-name="align_col"
#SBATCH -o logs/VCF_pipeline_%A_%a.out
#SBATCH -e logs/VCF_pipeline_%A_%a.err
#SBATCH --array=1-206
#SBATCH -N 1
#SBATCH -n 8
#"$SLURM_ARRAY_TASK_ID"

# activate our source bashrc and our conda environment
source /home/ssmith/.bashrc
source activate wgs_env

# get our sample name which corresponds with the SLURM_ARRAY_TASK_ID we are processing.
# This will pull the second column value (sample name) from the novo2name CSV file that matches with the SLURM_ARRAY_TASK_ID
sample_name=$(awk -v pattern="$SLURM_ARRAY_TASK_ID"  -F',' '$1 == pattern { print $2 }' /data2/ssmith/novo2name.csv)

# get the basenames of the sample we are processing (forward and reverse files)
file_name_1=$(basename /data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz)
file_name_2=$(basename /data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz)

echo BEGINNING FASTQC SCRIPT

# the directory path we need for fastQC
dirPath=/data2/ssmith/fastqs/

# run fastqc to carry out initial QC
fastqc "$dirPath""$SLURM_ARRAY_TASK_ID"_*.fq.gz --noextract -t 6 -a /data2/ssmith/adapters.txt --outdir=/data2/ssmith/QC/initial_QC/

echo TRIMMOMATIC ADAPTER TRIMMING

# run trimmomatic, see the doc for sepcific arguments
trimmomatic \
PE \
/data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz \
/data2/ssmith/fastqs/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz \
/data2/ssmith/trimmed/"$file_name_1" /data2/ssmith/discarded/"$file_name_1" \
/data2/ssmith/trimmed/"$file_name_2" /data2/ssmith/discarded/"$file_name_2" \
ILLUMINACLIP:/data2/ssmith/adapters.fa:2:30:10 \
MINLEN:50

echo BEGINNING FASTQC SCRIPT
# run fastQC on our trimmed reads
dirPath=/data2/ssmith/trimmed/
fastqc "$dirPath""$SLURM_ARRAY_TASK_ID"_*.fq.gz --noextract -t 6 -a /data2/ssmith/adapters.txt --outdir=/data2/ssmith/QC/after_QC/

echo BWA-MEM2 ALIGNMENT

# align the trimmed reads using bwa-mem2 and sort using samtools
bwa-mem2 mem -t 4 /data/ssmith/c_l_genome/apis_c_l_genome.fa  \
/data2/ssmith/trimmed/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_1.fq.gz \
/data2/ssmith/trimmed/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_2.fq.gz \
| samtools sort -@ 4 -o /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam

sortedPath=/data2/ssmith/bams

echo SAMTOOLS FLAGSTAT

# get mapping statistics using samtools
samtools flagstat "$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam \
> /data2/ssmith/QC/initial_QC/"$SLURM_ARRAY_TASK_ID"_"$sample_name".txt

echo PICARD MARK DUPLICATES 

# mark PCR duplicates
picard MarkDuplicates \
I="$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam \
M=/data2/ssmith/QC/after_QC/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_dup_metrics.txt \
O="$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.dupM.bam

rm "$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.bam

echo PICARD COLLEGE WGS METRICS

# collect mapping stats
picard CollectWgsMetrics \
I="$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.dupM.bam \
O=/data2/ssmith/QC/after_QC/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_collect_wgs_metrics.txt \
R=/data/ssmith/c_l_genome/apis_c_l_genome.fa

echo SAMTOOLS FLAGSTAT

# get new mapping stats from samtools after duplicate marking
samtools flagstat "$sortedPath"/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.dupM.bam \
> /data2/ssmith/QC/after_QC/"$SLURM_ARRAY_TASK_ID"_"$sample_name".txt

echo PICARD ADD OR REPLACE GROUPS

# add in sample ID using picard
picard AddOrReplaceReadGroups \
I=/data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.dupM.bam \
O=/data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".grpd.bam \
RGID=1 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=A_"$SLURM_ARRAY_TASK_ID"

rm /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".sorted.dupM.bam

# we cannot run gatk base recalibration as we do not have a known SNP file for the new reference genome
# indel recalibration is done through lofreq later

samtools index /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".grpd.bam

# indel quality control through lofreq
# the bam files created here are what we use later in the population genomics work (popoolation/poolfstat)
lofreq indelqual --dindel \
--ref /data/ssmith/c_l_genome/apis_c_l_genome.fa \
-o /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_indels.bam \
/data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".grpd.bam

rm /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".grpd.bam
rm /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name".grpd.bam.bai

samtools index /data2/ssmith/bams/"$SLURM_ARRAY_TASK_ID"_"$sample_name"_indels.bam

