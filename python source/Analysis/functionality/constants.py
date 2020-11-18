import numpy as np

DATA_ROOT = "D:\\SyntBio data\\response to alpha factor\\Roberts_2000_PMID_10657304\\"

SEQUENCE_MAPPING = {
    # DNA nucleotides
    'A': np.array([1.0, 0.0, 0.0, 0.0]),
    'C': np.array([0.0, 1.0, 0.0, 0.0]),
    'G': np.array([0.0, 0.0, 1.0, 0.0]),
    'T': np.array([0.0, 0.0, 0.0, 1.0]),

    # RNA nucleotide (bonus)
    'U': np.array([0.0, 0.0, 0.0, 1.0]),

    # 2-nucleotide combinations
    'W': np.array([0.5, 0.0, 0.0, 0.5]),
    'S': np.array([0.0, 0.5, 0.5, 0.0]),
    'M': np.array([0.5, 0.5, 0.0, 0.0]),
    'K': np.array([0.0, 0.0, 0.5, 0.5]),
    'R': np.array([0.5, 0.0, 0.5, 0.0]),
    'Y': np.array([0.0, 0.5, 0.0, 0.5]),

    # 3-nucleotide combinations
    'B': np.array([0.0, 1 / 3, 1 / 3, 1 / 3]),
    'D': np.array([1 / 3, 0.0, 1 / 3, 1 / 3]),
    'H': np.array([1 / 3, 1 / 3, 0.0, 1 / 3]),
    'V': np.array([1 / 3, 1 / 3, 1 / 3, 0.0]),

    # 4-nucleotide combinations
    'N': np.array([0.25, 0.25, 0.25, 0.25]),
}
