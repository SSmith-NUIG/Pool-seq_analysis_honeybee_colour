# Pool-seq_analysis_honeybee_colour
Pipeline to determine which regions/genes cause colour differences in honey bees  
In this analysis we have 3 honey bee colonies that display a mixed colour phenotype  
Half of the honey bees in the hive are pure black in colour and half have a yellow band on the abdomen  
We take two samples from each hive, 30 black and 30 banded for a total of 6 samples (3 hives * 2 samples per hive)  


## Step one
Run ```bam_creation.sh``` on the raw fastq files for all six samples  

## Step two
Run ```CMH_all_colour.sh``` to run CMH tests for our samples to find SNPs that are  
signficantly different between black and banded bees  

## Step three
Run ```cmh_results_analyis.py``` to combine the CMH results and create manhatten plots  
This also creates output files required for GOwinda analysis  

## Step four
Create the input file for GOwinda using ```gowinda_input_file.py```  
This final file was HEAVILY manually modified (days of work, contact me to obtain it)  
These manual modifications were necessary to ensure that the gene IDs in the input file matched the gene IDs in the GTF file.  

## Step five
Run ```gowinda_script.sh``` to obtain an output file that contain information for each signficiant metabolic process 
discovered. This can then be searched for processess involved in pigmentation/development/melanin etc.



