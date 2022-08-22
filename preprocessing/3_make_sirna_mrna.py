import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from params

### create a sirna rev-compl --> mRNA target file

records = list(SeqIO.parse(params.mrna_fasta_file, "fasta"))
sirna = params.sirna_fasta_file
array = []
count = 0

for seq_record in SeqIO.parse(sirna,'fasta'):
# for sirin, seqin in sirna.iterrows():
    count += 1
    sir = seq_record.seq.upper().replace('U', 'T')
    sir_id = seq_record.id
    # sir = seqin[1].upper().replace('U', 'T')
    # sir = Seq(sir)
    rev = Seq.reverse_complement(sir)
    for r in records:
        if rev in r.seq:
            array.append([sir_id, rev, r.id, r.seq])


pd.DataFrame(array).to_csv('datase_tofold.csv', index=False,header = False)
