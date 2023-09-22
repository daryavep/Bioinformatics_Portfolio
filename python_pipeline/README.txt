#RBIF100 Week 6 Assignment
#Author: Darya Veprinski

The goal of this script is to analyze mold data and identify the mutations that resulted in the different color molds.

pipeline.py is divided into 5 steps:
1. Demultiplexing the pooled fastq and trimming the end of the lower quality reads.
2. Permorming alignment on each FASTQ to the reference sequence and generating samfiles.
3. Converting the samfiles to bamfiles.
4. Discovering SNPs.
5. Generating a report.

###########################################################################################
Execute the script inside of the week6 folder by typing the following command:

python3 pipeline.py --fastq hawkins_pooled_sequences.fastq

###########################################################################################

Please allow a few seconds for the script to iterate through all the files. When it is done, it will print:
"All done! Please review the following contents: directory 'fastqs', directory 'sorted_bam', and file 'report.txt'."


The necessary files for the script are located in week6 under the names:
dgorgon_reference.fa
harrington_clinical_data.txt
hawkins_pooled_sequences.fastq

This script will create:
1. directory fastqs: contains individual fastq files for each sample.
2. directory sorted_bam: containes the sorted and indexed bam files.
3. file report.txt: this is the final report for the findings.
4. directory bwa_index: keeps bwa products organized, can be disregarded.
5. file mutations-report.txt: keeps mutation reports to be populated into step five's report, can be disregarded.