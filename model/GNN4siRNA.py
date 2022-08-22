#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 13:10:50 2022

@author: fiannaca
"""



import stellargraph as StellarGraph
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
import scipy
import scipy.sparse
import scipy.sparse.linalg

from stellargraph.mapper import HinSAGENodeGenerator
from stellargraph.layer import HinSAGE
from tensorflow.keras import layers, Model, optimizers


# import file with parameters
import params



# Specify the minibatch size and the number of epochs for training the model
batch_size = params.batch_size
epochs = params.epochs


#############################################
# import file with sirna / target thermodynamic features
#############################################

# k-mers of siRNA sequences
sirna_pd = pd.read_csv( params.sirna_kmer_file, header=None )
sirna_pd= sirna_pd.set_index(0)

# k-mers of mRNA sequences
mRNA_pd = pd.read_csv( params.mrna_kmer_file, header=None )
mRNA_pd= mRNA_pd.set_index(0)

# thermodynamic features of siRNA-mRNA interaction
thermo_feats_pd = pd.read_csv( params.sirna_target_thermo_file, header=None )

# sirna_efficacy_values
sirna_efficacy_pd = pd.read_csv( params.sirna_efficacy_file )


# rename first 2 columns in "source" and "target"
thermo_feats_pd.rename(columns={0 : "source", 1 : "target"}, inplace= True)

# Here we transform interaction edges in "interaction nodes"
# Intercation node has 2 edges that connect it to siRNA and mRNA, respectively
# Node ID cames from source and target ids
interaction_pd = thermo_feats_pd.drop(['source','target'], axis=1)
interaction_pd["index"] = thermo_feats_pd['source'].astype(str) + "_" + thermo_feats_pd['target']
interaction_pd= interaction_pd.set_index("index")

# New edges have no features
sirna_edge_pd_no_feats = thermo_feats_pd[['source','target']]
data1 = { 'source' : list(interaction_pd.index),
         'target' : sirna_edge_pd_no_feats['source']}
data2 = { 'source' : list(interaction_pd.index),
         'target' : sirna_edge_pd_no_feats['target']}

all_my_edges = pd.DataFrame(data1)
all_my_edges_temp = pd.DataFrame(data2)

# Merge all the edges 
all_my_edges = pd.concat([all_my_edges, all_my_edges_temp], ignore_index = True, axis = 0)


# We want to predict the interaction weight, i.e. the label of interaction node
interaction_weight = sirna_efficacy_pd['efficacy']
interaction_weight = interaction_weight.set_axis(interaction_pd.index)



# Create Stellargraph object
my_stellar_graph = StellarGraph.StellarGraph( {"siRNA": sirna_pd, "mRNA": mRNA_pd, "interaction": interaction_pd}, 
                                             edges=all_my_edges, source_column="source", target_column= "target") 



################################################
# Create the model
################################################

overall_PCC = []
overall_mse = []

# HinSAGE parameters
hinsage_layer_sizes = params.hinsage_layer_sizes
hop_samples = params.hop_samples
dropout = params.dropout
loss_function= params.loss
learning_rate= params.lr



# with range = 1, it make only a repeat
for repeat in range(1):
    score_PCC = []
    score_mse= []
    myfold = 0
    
    # perform 10-fold cross validation
    kfold = KFold(n_splits=10, shuffle=True, random_state = 2)
    for train, test in kfold.split(interaction_weight):
        train_interaction, test_interaction = interaction_weight[train], interaction_weight[test]

        myfold = myfold + 1
    
        print("Repeated n. ", repeat, "-- Fold n. ", myfold)


# Create the generators to feed data from the graph to the Keras model
# We specify we want to make node regression on the "interaction" node
        generator = HinSAGENodeGenerator(
            my_stellar_graph, batch_size, hop_samples, head_node_type= "interaction")

        train_gen = generator.flow(train_interaction.index, train_interaction, shuffle=True)


        hinsage_model = HinSAGE(
            layer_sizes=hinsage_layer_sizes, generator=generator, bias=True, dropout=dropout
            )

# Expose input and output sockets of hinsage:
        x_inp, x_out = hinsage_model.in_out_tensors()

        prediction = layers.Dense(units=1)(x_out)


# Now let’s create the actual Keras model with the graph inputs x_inp 
# provided by the graph_model and outputs being the predictions
        model = Model(inputs=x_inp, outputs=prediction)
        model.compile(
            optimizer=optimizers.Adam(lr=learning_rate),
            loss=loss_function
            )


# Train the model, keeping track of its loss and accuracy on the training set, 
# and its generalisation performance on the test set

        test_gen = generator.flow(test_interaction.index, test_interaction)


        history = model.fit(
            train_gen, epochs=epochs, validation_data=test_gen, verbose=2, shuffle=False
            )


# Plot the training history, epochs Vs loss
#        StellarGraph.utils.plot_history(history)


# Now we have trained the model we can evaluate on the test set.
        pred= model.predict(test_gen)

# Calculate the PCC and MSE scores
        pearson = scipy.stats.pearsonr(test_interaction,pred)
        score_PCC.append(pearson[0])
    
        mse_run = mean_squared_error(test_interaction, pred)
        score_mse.append(mse_run)
        
    

    overall_PCC.append(np.mean(score_PCC))
    overall_mse.append(np.mean(score_mse))


# Print the average values on PCC and MSE
print ("overall PCC score = ", np.mean(overall_PCC))
print ("overall MSE score = ", np.mean(overall_mse))