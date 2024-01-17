#!/bin/sh

java -Xmx2g -jar /home/stephen/gowinda/Gowinda-1.12.jar \
--snp-file /home/stephen/gowinda/fixing_input/triple_test/all_snps_cmh.csv \
--candidate-snp-file /home/stephen/gowinda/fixing_input/triple_test/gowinda_significant_all_genome_fixed_new.csv \
--annotation-file /home/stephen/gowinda/fixing_input/triple_test/gowinda_gtf3.gtf \
--gene-set-file /home/stephen/gowinda/fixing_input/triple_test/gowinda_add_HGD2.csv \
--output-file /home/stephen/gowinda/fixing_input/triple_test/gowinda_snp_results_all_genometesting_fixed.txt \
--simulations 50000 --min-significance 1 \
--gene-definition updownstream2000 \
--mode snp \
--detailed-log \
--threads 32 \
--min-genes 1
