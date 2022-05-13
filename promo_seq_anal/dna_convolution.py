import numpy as np
from Bio.Seq import Seq
from promo_seq_anal.str_to_onehot import dna_to_onehot

def conv_1d(seq, mat, dna_len, kernal_size):
    res = []
    for i in range(dna_len):
        dna_cut = seq[:, i:i + kernal_size]
        conv = np.sum(dna_cut * mat, axis=0)
        conv = np.prod(conv)
        res.append(conv)

    return res


def get_dna_convolution(dna, mat):
    dna_len = len(dna)
    dne_rev_compl = str(Seq(dna).reverse_complement())

    dna = dna_to_onehot(dna)
    dne_rev_compl = dna_to_onehot(dne_rev_compl)

    if mat.shape[1] % 2 == 0:
        mat = np.concatenate([mat, np.ones((4, 1), 'float64')], dtype='float64', axis=-1)

    ovelap_size = mat.shape[1] // 2
    kernal_size = mat.shape[1]

    dna = np.concatenate([np.zeros((4, ovelap_size)), dna, np.zeros((4, ovelap_size))], dtype='float64', axis=-1)
    dne_rev_compl = np.concatenate([np.zeros((4, ovelap_size)), dne_rev_compl, np.zeros((4, ovelap_size))],
                                   dtype='float64', axis=-1)

    result_ = conv_1d(dna, mat, dna_len, kernal_size)
    result_complement = conv_1d(dne_rev_compl, mat, dna_len, kernal_size)

    result_ = np.array(result_)[None, :]
    result_complement = np.flip(np.array(result_complement))[None, :] * -1

    result = np.concatenate([result_, result_complement], axis=0)

    return result


if __name__ == '__main__':
    pCLN1 = 'GAATGTTACGGGACTAACAGGCGGATGTAAATTACTCACTTAAACGCAGCCAAACATCATCGAGAACTTAGGGTAGCGTGCCACAAAATTTGCATGAATAAACTTTTGTTTTCCTAATTCGACAGCATTCCCTTGTTCGCAACACTTCACTGATAGGAAATCGAATAGCGCACACTCTCTTCTGGGACATACCCCAATGCGGTAAAGCCACGAAAACACCGCGCGTAAAGGGGTAAACAAGTCCATTCCTACAACCTCTTGGAGAAATTCTTTACCTACTACAACCCCGCGCCTGATACTTTCAGTATTCATGACAACTCGAGCCAGATCCCGCTCGTGGGCGTGTTCATTCTGTGACGATCCACTAGCGACTTCTTTGTTCAGCCTGCAAGAGACGCGTTCAAGGAAGAATTCGCGATTTTACTTCTTCGAGGGAATCTCGCACCGCGTTAGTTAGTTTCCAACCTTGAAAGCATCGGAGACGCATTTTTGGCGATTTTGCTGGATTGAGCTGAATGGTGCCAGGTCGAGGCTGGGAGGGAGACTAACTCGAAAGTGACGAAGACTCGAAAATTAAAAAAAAAGATACTGCAGAAGGCAAGATTGAGAATGGAGTAAAGGCAGCGTGGGTCCCCTGTGGAAACCGCAGTTTTCCTGCGCCAAGTGGTACCGGTGCGAGTGCAGCAATTAATCTCTCGATATTTTCTTAGTATCTCTTTTTATATAAGAATATATTTTGGAATTGGTAATGCTTATCTTCAATAGTTTCTTAGTTGAATGCACACTTAAGAGCAAATTGGCCAAGGAGTTCTTCGTTCGCTTTAATTTATTTCCTGGTTATTGTCAATTTATTCATCCCATCTCCCCAGGATAGAAGAAATTAGTGTAATTTTGCTGACAATACATTTTAACGACGATAACAATAATAGCAATTAAATAAAATAGCACTACCACCACTCCACTGCTCGTTAGCTATTTCTGTAAAATAAATAAAAAGATC'



