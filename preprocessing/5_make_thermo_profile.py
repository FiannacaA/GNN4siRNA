# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:53:19 2022

@author: Massimo La Rosa
crate the complete thermodynamic profile file, merging thermodynamic features 
and folding free energy
"""

import pandas as pd
import params

path_thermo = "dataset_th.csv" #thermodynamic features
path_sirna = params.sirna_mrna_efficacy #sirna-mrna-efficacy file
path_rnaup = "dataset_folded.csv" #folding features

thermo = pd.read_csv(path_thermo, header=0)
sirna = pd.read_csv(path_sirna, header=0)
rnaup = pd.read_csv(path_rnaup, header=0)

df = pd.DataFrame()
df = pd.concat([df,sirna.iloc[:,0:2]],axis=1,)
df = pd.concat([df,thermo.iloc[:,1:20]],axis=1,)
df = pd.concat([df,rnaup.iloc[:,2:5]],axis=1)

df.to_csv("sirna_target_thermo.csv",header = False,index=False)
                