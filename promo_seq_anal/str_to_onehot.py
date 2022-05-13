import numpy as np

NUCLEOTIDE_MAP = {
    'A': 0,
    'T': 1,
    'G': 2,
    'C': 3
}

NUCLEOTIDE_MAP_ = {
    0: 'A',
    1: 'T',
    2: 'G',
    3: 'C'
}


def dna_to_onehot(dna):
    seq_digits = [NUCLEOTIDE_MAP[nuc] for nuc in dna]

    seq_onehot = np.zeros((4, len(dna)), 'float64')
    seq_onehot[seq_digits, np.arange(len(dna))] = 1.0

    return seq_onehot
