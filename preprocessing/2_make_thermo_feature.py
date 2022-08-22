# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:53:19 2022

@author: Massimo La Rosa
compute the thermodynamic features
"""
import pandas as pd
from Bio import SeqIO
import params


def build_kmers(sequence, ksize):
    kmers = []
    n_kmers = len(sequence) - 1

    for i in range(n_kmers):
        kmer = sequence[i:i + ksize]
        kmers.append(kmer)
    return kmers



path_sirna = params.sirna_fasta_file #input fasta file
total_compute = []
intermolecular_initiation = 4.09
simmetry_correction = 0.43

# calculation of stability based on bimers on the sequence
# for index, sequence in df.iterrows():
for seq_record in SeqIO.parse(path_sirna,'fasta'):
    sum = 0
    single_sum = []
    
    sequence = seq_record.seq
    sequence = sequence[:19]
    seq_id = seq_record.id
    bimer = build_kmers(sequence, 2)
    single_sum.append(seq_id)

    # checking if first is 'A' and last is 'U'
    if sequence[0] == 'A':
        sum += +0.45
    if sequence[18] == 'U':
        sum += +0.45

    # checking the kind of bimers
    for b in bimer:
        if b == 'AA' or b == 'UU':
            single_sum.append(-0.93)
            sum+=-0.93
        elif b == 'AU':
            single_sum.append(-1.10)
            sum += -1.10
        elif b == 'UA':
            single_sum.append(-1.33)
            sum += -1.33
        elif b == 'CU' or b == 'AG':
            single_sum.append(-2.08)
            sum += -2.08
        elif b == 'CA' or b == 'UG':
            single_sum.append(-2.11)
            sum += -2.11
        elif b == 'GU' or b == 'AC':
            single_sum.append(-2.24)
            sum += -2.24
        elif b == 'GA' or b == 'UC':
            single_sum.append(-2.35)
            sum += -2.35
        elif b == 'CG':
            single_sum.append(-2.36)
            sum += -2.36
        elif b == 'GG' or b == 'CC':
            single_sum.append(-3.26)
            sum += -3.26
        elif b == 'GC':
            single_sum.append(-3.42)
            sum += -3.42
        else:
            single_sum.append(0)
            sum += 0

    sum += intermolecular_initiation
    sum += simmetry_correction
    single_sum.append(sum)
    total_compute.append(single_sum)

pd.DataFrame(total_compute).to_csv('dataset_th.csv', index=False, doublequote=False)
