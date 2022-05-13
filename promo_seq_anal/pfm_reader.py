import pandas as pd
import numpy as np


def read_pfm(f_path, ignore_sum=False):
    data = pd.read_csv(f_path, delimiter='\t', index_col=False, header=None)

    assert data.shape[0] == 4 and data.shape[1] > 0 and ''.join(list(data[0])) == 'ATGC'

    pfm = np.array(data.loc[:, 1:])

    assert ignore_sum or np.sum(np.abs(np.sum(pfm, axis=0) - np.ones((pfm.shape[1])))) < (1e-3 * pfm.shape[1])

    return pfm
