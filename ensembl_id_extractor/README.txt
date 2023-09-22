#RBIF100 Week 8 Assignment
#Author: Darya Veprinski

The goal of this script is to extract any gene ID from Ensembl using RESTful API, and then generate homologies.

###########################################################################################
Execute the script inside of the week8 folder by typing the following command:

python3 gene_homologs.py MC1R

###########################################################################################

NOTE: You may insert any gene name in place of MC1R for different results.

There are no additional files necessary for executing the script.

This script will create:
1. A printed statement of the selected gene Ensembl ID. 
2. A FASTA file containing the gene's DNA sequence and amino acid sequence.
3. A .txt file containing the gene's unique homologs.