# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 10:53:19 2022

@author: Massimo La Rosa
compute the thermodynamic folding free energy by means of RNAUp tool of
ViennaRNA package - https://www.tbi.univie.ac.at/RNA/
"""
import subprocess
import pandas as pd

# python bash script to calculate all the interactions
stability = pd.read_csv("datase_tofold.csv", header=None)
array = []

for ind, seq in stability.iterrows():
    couple = seq[3] + '\n' + seq[1]
    # calling RNAup with some specifics
    proc = subprocess.Popen(["RNAup", "-b", "-d2", "--noLP", "-c", "'S'", "RNAup.out"],
                            stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        outs, errs = proc.communicate(input=couple.encode(encoding='utf-8'))
        start = ' ('
        end = ')\n'
        s = outs.decode(encoding='utf-8')
        s = s[s.find(start) + len(start):s.rfind(end)]
        s = s.replace('=', '').replace('+', '')
        s = list(map(float, s.split()))
        array.append([seq[0], s[0], s[2], s[3]])

    except subprocess.TimeoutExpired:
        proc.kill()
        outs, errs = proc.communicate()

# printing the result into csv
pd.DataFrame(array).to_csv('dataset_folded.csv')
