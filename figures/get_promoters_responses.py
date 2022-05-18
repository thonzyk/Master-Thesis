import numpy as np
from matplotlib import pyplot as plt
from Bio.Seq import Seq
import json

import pandas as pd
import cufflinks as cf
import chart_studio.plotly as py
import seaborn as sns
import plotly.express as px

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


def get_sequence_response(seq, kernel):
    assert len(kernel.shape) == 2
    assert seq.shape[0] == kernel.shape[0]

    seq_len = seq.shape[1]
    kernel_size = kernel.shape[1]

    response = []

    seq_dec = decode_sequence(seq)

    for i in range(seq_len - kernel_size + 1):
        seq_part = seq[:, i:i + kernel_size]

        this_resp = np.product(np.sum(seq_part * kernel, axis=0))
        if this_resp > 0.25:
            seeec = seq_dec[i:i + kernel_size]
            print(seeec, '\t', Seq(seeec).reverse_complement(), '\t', this_resp, '\t', i)
        response.append(this_resp)

    return np.array(response)


def get_tec1_promoters_responses():
    with open('data/tec1_recognition.tsv', 'r', encoding='utf-8') as fr:
        rows = []
        for line in fr:
            items = line[:-1].split('\t')
            items = np.array([float(item) for item in items])[None, :]
            rows.append(items)

        tec1_mat = np.concatenate(rows, axis=0)

    with open('data/CLN1_F.txt', 'r', encoding='utf-8') as fr:
        cln1_f = 'TCTGGA'
    with open('data/CLN1_R.txt', 'r', encoding='utf-8') as fr:
        cln1_r = fr.readline()[:-1]
    with open('data/GAS1_F.txt', 'r', encoding='utf-8') as fr:
        gas1_f = fr.readline()[:-1]
    with open('data/GAS1_R.txt', 'r', encoding='utf-8') as fr:
        gas1_r = fr.readline()[:-1]

    cln1_f, cln1_r, gas1_f, gas1_r = encode_sequence(cln1_f), encode_sequence(cln1_r), encode_sequence(
        gas1_f), encode_sequence(gas1_r)

    print('=== CLN1-F ===')
    cln1_f_resp = get_sequence_response(cln1_f, tec1_mat)
    print('=== CLN1-R ===')
    cln1_r_resp = -get_sequence_response(cln1_r, tec1_mat)
    print('=== GAS1-F ===')
    gas1_f_resp = get_sequence_response(gas1_f, tec1_mat)
    print('=== GAS1-R ===')
    gas1_r_resp = -get_sequence_response(gas1_r, tec1_mat)

    cln1_resp = list(np.concatenate([cln1_f_resp, cln1_r_resp]))
    gas1_resp = list(np.concatenate([gas1_f_resp, gas1_r_resp]))

    resp_len = len(cln1_f_resp)
    f_x_tics = [*range(0, resp_len)]
    r_x_tics = [*range(resp_len, 0, -1)]
    x_tics = f_x_tics + r_x_tics

    cln1 = {'x': x_tics, 'y': cln1_resp}
    gas1 = {'x': x_tics, 'y': gas1_resp}

    with open('data/cln1.json', 'w', encoding='utf-8') as fw:
        json.dump(cln1, fw, indent=4)

    with open('data/gas1.json', 'w', encoding='utf-8') as fw:
        json.dump(gas1, fw, indent=4)

    plt.figure()
    plt.plot(x_tics, cln1_resp)

    plt.figure()
    plt.plot(x_tics, gas1_resp)

    plt.show()


def plotly_responses():
    with open('data/cln1.json', 'r', encoding='utf-8') as fr:
        cln1 = json.load(fr)
    with open('data/gas1.json', 'r', encoding='utf-8') as fr:
        gas1 = json.load(fr)

    arr_1 = np.random.rand(50, 4)
    df_1 = pd.DataFrame(arr_1, columns=['A', 'B', 'C', 'D'])
    df_1.iplot()


if __name__ == '__main__':
    get_tec1_promoters_responses()
