# msa_conservation_index
Per-position conservation index of multiple sequence alignment

### How to run
```
msa_conservation_index.py msa.fa
```

### Input and output formats

Input data is a multiple sequence alignment (MSA) in fasta format.

Output goes to standard out, with 3 tab-delimited columns:
- 1-based position (column) within the MSA.
- A string with the amino acids found in such position.
- A conservation index (Ci) computed for such position.

#### Output example

```
...
234     SSE-E--DDDDDDEEEEEEE    0.5425222162502829
235     LLL-L--LLLLLLLSSSSSS    0.4750000000000001
236     PPP-Q--PPPPPPLLQTRRR    0.3652547316917079
237     NNN-N--NNNNNNNNNNNNN    0.7224999999999999
238     EEK-K--KKKRRRRRRRRRR    0.5613888888888889
239     KKK-V--KKKKKKKKKKKKK    0.6613561808316413
240     TTT-T--TTTTTTTTTTTTT    0.7224999999999999
241     VVI-L--VVVVVVVVVVVVV    0.6837499999999999
242     KKR-R--RRRRRRRRRRRRR    0.6725
243     III-I--VVVVVVIIIVVIV    0.6775000000000001
...
```

### Protein positions

The command outputs by default MSA positions (see above).
However, positions for each protein are different, due to the presence of gaps.
To add additional columns showing the (1-based) position of each Ci for each protein, run the command using a second argument.
For example:
```
msa_conservation_index.py msa.fa true
```
The new columns will correspond to the sequences in the same order as they appear in the MSA.
For example, the 4th column will show the positions for the 1st protein in the MSA, the 5th column will show the positions for the 2nd protein, and so on.

### How it works

Based on https://doi.org/10.1093/bioinformatics/17.8.700

- load blosum62 substitution matrix and shift it from negative-positive to only non-negative values                                           
- normalize amino acids scores of the non-negative blosum62 with Sab=Sab/squareroot(Saa*Sbb)                                              
- At each MSA position:                                                                            
-- Compute unweighted aas frequencies                                                             
-- Compute "sum-of-pairs" conservation index (Ci)                                                 
- Compute per-position Ci for each protein.
- Normalize Ci to (0,1) by dividing by the number of rows (sequences), taking into account gaps.
- Output per-position normalized Ci.

