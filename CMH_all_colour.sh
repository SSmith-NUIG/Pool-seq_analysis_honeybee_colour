#!/bin/sh 
#SBATCH --job-name="CMH_cell"
#SBATCH -o /data/ssmith/logs/CMH_cell_%A_%a.out
#SBATCH -e /data/ssmith/logs/CMH_cell_%A_%a.err
#SBATCH -N 1
#SBATCH -n 4
#"$SLURM_ARRAY_TASK_ID"
# 33-38,40-57%12

#source activate wgs_env


#/home/ssmith/samtools-0.1.16/samtools mpileup -f /data/ssmith/c_l_genome/apis_c_l_genome.fa -q 20 -Q 20 -B \
#-b /data/ssmith/scripts/colony_analysis/jan_analysis/pooled_pop_pipeline/colour_analysis/colours_bamlist.txt \
#-l /data/ssmith/scripts/colony_analysis/jan_analysis/pooled_pop_pipeline/cl.bed \
#> /data3/ssmith/ena/population_genomics/colour_all.mpileup

#module load java

#java -jar /home/ssmith/popoolation2_1201/mpileup2sync.jar \
#--input /data3/ssmith/ena/population_genomics/colour_all.mpileup \
#--output /data3/ssmith/ena/population_genomics/colour_all.sync \
#--min-qual 20 \
#--fastq-type sanger \
#--threads 12

module load R/R-4.0.2

perl /home/ssmith/popoolation2_1201/cmh-test.pl \
--input /data3/ssmith/ena/population_genomics/colour_all.sync \
--output /data3/ssmith/ena/population_genomics/colour_1_1_2_2_3_3.cmh \
--min-count 10 \
--min-coverage 15 \
--max-coverage 200 \
--population 1-2,3-4,5-6

perl /home/ssmith/popoolation2_1201/cmh-test.pl \
--input /data3/ssmith/ena/population_genomics/colour_all.sync \
--output /data3/ssmith/ena/population_genomics/colour_1_2_2_3_3_1.cmh \
--min-count 10 \
--min-coverage 15 \
--max-coverage 200 \
--population 1-4,3-6,5-2

perl /home/ssmith/popoolation2_1201/cmh-test.pl \
--input /data3/ssmith/ena/population_genomics/colour_all.sync \
--output /data3/ssmith/ena/population_genomics/colour_1_3_2_1_3_2.cmh \
--min-count 10 \
--min-coverage 15 \
--max-coverage 200 \
--population 1-6,3-2,5-4
