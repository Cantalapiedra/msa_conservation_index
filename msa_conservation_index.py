#!/usr/bin/env python3
## CPCantalapiedra 2021

import math, sys
from collections import Counter

from Bio.Align import substitution_matrices
from Bio import AlignIO

##
# Arguments

msafn = sys.argv[1]

##
# Functions

def load_substitution_matrix():
    
    ##
    # Prepare the substitution matrix

    # Load matrix

    # names = substitution_matrices.load()
    # print(names)
    matrix = substitution_matrices.load('BLOSUM62')

    # Shift matrix to minimum value = 0

    matrix_min = min(matrix.values())
    matrix = matrix - matrix_min

    # Normalize blosum62 with Sab = Sab/squareroot(Saa*Sbb)
    # according to https://doi.org/10.1093/bioinformatics/17.8.700

    norm_mat = matrix.copy()
    for (aa1, aa2),Sab in matrix.items():
        Saa = matrix[(aa1, aa1)]
        Sbb = matrix[(aa2, aa2)]
        Sab_2 = Sab / math.sqrt(Saa * Sbb)
        norm_mat[(aa1, aa2)] = Sab_2

    return norm_mat

##
# Main

# Prepare substitution matrix

norm_mat = load_substitution_matrix()

# Load MSA

msa = AlignIO.read(msafn, "fasta")

# Compute per-position conservation index (Ci)

msa_num_seqs = len(msa)
msa_num_cols = msa.get_alignment_length()

for col in range(msa_num_cols):
    col_aas = msa[:, col]
    
    # Compute unweighted aas frequencies
    # according to https://doi.org/10.1093/bioinformatics/17.8.700

    col_counter = Counter(col_aas)
    freqs = {aa:(col_counter[aa] / msa_num_seqs) for aa in col_counter if aa != "-"}

    # Compute "sum-of-pairs" conservation index (Ci)
    # according to https://doi.org/10.1093/bioinformatics/17.8.700

    Ci = 0
    for aa1,f1 in freqs.items():
        for aa2,f2 in freqs.items():
            Ci += f1*f2*norm_mat[(aa1, aa2)]
            
    print(f"{col}\t{col_aas}\t{Ci}")

        

## END
