import numpy as np
from Bio.Seq import Seq

ENCODER = {
    'A': 0,
    'T': 1,
    'G': 2,
    'C': 3
}

DECODER = {ENCODER[key]: key for key in ENCODER}


def encode_sequence(seq):
    seq = seq.upper()
    encoded_seq = []
    for letter in seq:
        assert letter in ENCODER
        apendix = np.zeros((4, 1), 'float64')
        apendix[ENCODER[letter], 0] = 1.0
        encoded_seq.append(apendix)

    encoded_seq = np.concatenate(encoded_seq, axis=1)
    return encoded_seq


def decode_sequence(seq):
    seq = np.argmax(seq, axis=0)
    seq_decoded = ""
    for token in list(seq):
        seq_decoded += DECODER[token]

    return seq_decoded


def get_score(seq, rbam):
    seq = encode_sequence(seq)
    resp = seq * rbam
    resp = np.sum(resp, axis=0)
    resp = np.prod(resp)

    return resp


def get_seq_response(seq_f, rbam):
    seq_len = len(seq_f)
    rbam_len = rbam.shape[1]

    score_99 = 0
    score_50 = 0
    score_20 = 0

    response_f = []
    response_r = []
    for i in range(seq_len - rbam_len):
        subseq_f = seq_f[i:i + rbam_len]
        subseq_r = str(Seq(subseq_f).reverse_complement())

        score_f = get_score(subseq_f, rbam)
        score_r = get_score(subseq_r, rbam)

        response_f.append(score_f)
        response_r.append(score_r)

        score_99 += int(score_f > 0.98) + int(score_r > 0.98)
        score_50 += int(score_f > 0.49) + int(score_r > 0.49)
        score_20 += int(score_f > 0.19) + int(score_r > 0.19)

    assert score_99 <= score_50 <= score_20

    return {
        'response_f': response_f,
        'response_r': response_r,
        'score_99': score_99,
        'score_50': score_50,
        'score_20': score_20
    }


def load_binding_mat(path='source_data/Heise_2010_PMID_20118212/tec1_recognition.tsv'):
    mat = []
    with open(path, 'r', encoding='utf-8') as fr:
        for line in fr:
            items = line[:-1].split('\t')
            items = [float(el) for el in items[1:]]
            mat.append(np.array(items)[None, :])

    mat = np.concatenate(mat, axis=0)
    return mat
