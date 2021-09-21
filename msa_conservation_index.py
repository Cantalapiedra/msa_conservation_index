#!/usr/bin/env python3
## CPCantalapiedra 2021

import math, sys
from collections import Counter

from Bio.Align import substitution_matrices
from Bio import AlignIO

##
# Arguments

msafn = sys.argv[1]
if len(sys.argv) > 2:
    show_per_prot_position = True
else:
    show_per_prot_position = False

##
# Functions

#
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

#
def compute_freqs(col_aas, msa_num_seqs):
    '''
    ComputeS unweighted aas frequencies
    according to https://doi.org/10.1093/bioinformatics/17.8.700
    '''
    freqs = None
    
    col_counter = Counter(col_aas)
    freqs = {aa:(col_counter[aa] / msa_num_seqs) for aa in col_counter if aa != "-"}
    
    return freqs

#
def compute_Ci(freqs, norm_mat):
    '''
    Compute "sum-of-pairs" conservation index (Ci)
    according to https://doi.org/10.1093/bioinformatics/17.8.700
    '''
    Ci = 0
    
    for aa1,f1 in freqs.items():
        for aa2,f2 in freqs.items():
            Ci += f1*f2*norm_mat[(aa1, aa2)]

    return Ci

#
def update_per_prot_positions(col_aas, prots_positions):
    prots_positions_str = []
    for i,aa in enumerate(col_aas):
        if aa != "-":
            prots_positions[i] += 1
            prots_positions_str.append(str(prots_positions[i]))
        else:
            prots_positions_str.append("-")
            
    return "\t".join(prots_positions_str)


##
# Main

# Prepare substitution matrix

norm_mat = load_substitution_matrix()

# Load MSA

msa = AlignIO.read(msafn, "fasta")

msa_num_seqs = len(msa)
msa_num_cols = msa.get_alignment_length()

# Prepare a dict to track position of each protein
# avoiding gaps. The key is the number of sequence
# in the alignment (1st sequence has key 1,
# 2nd sequence has key 2, etc...)

prots_positions = {i:0 for i in range(msa_num_seqs)}

# Prepare and output header

if show_per_prot_position == True:
    print("msa_pos\talignment\tCi", end = "")
    for i in range(0, msa_num_seqs - 1):
        print(f"\t{msa[i,:].id}", end = "")
    print()
else:
    print("msa_pos\talignment\tCi")


# Compute per-position conservation index (Ci)

for col in range(msa_num_cols):
    col_aas = msa[:, col]

    # Compute freqs

    freqs = compute_freqs(col_aas, msa_num_seqs)

    # Compute conservation index (Ci)

    Ci = compute_Ci(freqs, norm_mat)
    
    # Update position for each protein

    if show_per_prot_position == True:
        prots_positions_str = update_per_prot_positions(col_aas, prots_positions)
        print(f"{col+1}\t{col_aas}\t{Ci}\t{prots_positions_str}")
        
    else:
        print(f"{col+1}\t{col_aas}\t{Ci}")

## END
