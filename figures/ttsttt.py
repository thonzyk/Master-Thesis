import numpy as np

if __name__ == '__main__':
    with open('data/tec1_recognition.tsv', 'r', encoding='utf-8') as fr:
        rows = []
        for line in fr:
            items = line[:-1].split('\t')
            items = np.array([float(item) for item in items])[None, :]
            rows.append(items)

        tec1_mat = np.concatenate(rows, axis=0)

    S = np.array([
        [0, 0, 0, 0, 0, 1],
        [1, 0, 1, 0, 0, 0],
        [0, 0, 0, 1, 1, 0],
        [0, 1, 0, 0, 0, 0]
    ])

    res = np.matmul(np.transpose(S), tec1_mat)
    res = np.diag(res)
    res = np.prod(res)

    print()
