##################
# INPUT FILES: path of the two fasta files
##################

# siRNA sequences
sirna_fasta_file = "../data/raw/dataset_1/sirna_1.fas"
# mRNA sequences
mrna_fasta_file = "../data/raw/dataset_1/mRNA_1.fas"
# siRNA-mRNA efficacy list
sirna_mrna_efficacy = "../data/raw/dataset_1/sirna_mrna_efficacy.csv"

##################
# Preprocessing parameters
##################

# k value for k-mer representation
k_sirna = 3
k_mrna = 4