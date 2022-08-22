##################
# INPUT FILES: path of the three input files
##################

#k-mers of siRNA sequences
sirna_kmer_file = "../data/processed/dataset_2/sirna_kmers.txt"
#k-mers of mRNA sequences
mrna_kmer_file = "../data/processed/dataset_2/target_kmers.txt"
#thermodynamic features of siRNA-mRNA interaction
sirna_target_thermo_file = "../data/processed/dataset_2/sirna_target_thermo.csv"
# sirna_efficacy_values
sirna_efficacy_file = "../data/processed/dataset_2/sirna_mrna_efficacy.csv"

##################
# GNN Hyperparameters
##################

# specify the minibatch size and the number of epochs for training the GNN model
batch_size = 60
epochs = 10

# two hidden layers HinSAGE sizes
hinsage_layer_sizes = [32, 16]

# sizes of 1- and 2-hop neighbour samples for each hidden layer of the HinSAGE model
hop_samples = [8, 4]

# Dropout value for the HinSAGE model
dropout = 0.15

# Adamax learning rate
lr= 0.001

# loss function
loss='mse'




