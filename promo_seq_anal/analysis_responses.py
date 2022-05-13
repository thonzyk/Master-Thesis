import pandas as pd
from matplotlib import pyplot as plt
import os

from promo_seq_anal.constants import *

files = os.listdir(RESULT_DIR / 'promoter_tf_response' / 'TEC1')

for i, f_name in enumerate(files):
    print(i, ' ', f_name)
    if not f_name.endswith('.tsv'):
        continue
    data_tec = pd.read_csv(RESULT_DIR / 'promoter_tf_response' / 'TEC1' / f_name, delimiter='\t', index_col=False)
    data_ste = pd.read_csv(RESULT_DIR / 'promoter_tf_response' / 'STE12_consensus' / f_name, delimiter='\t',
                           index_col=False)

    plt.plot([*range(len(data_tec))], data_tec['response'], c='tab:blue')
    plt.plot([*range(len(data_tec))], data_tec['response_rev'], c='tab:blue')

    plt.plot([*range(len(data_ste))], data_ste['response'], c='tab:orange')
    plt.plot([*range(len(data_ste))], data_ste['response_rev'], c='tab:orange')

    plt.legend(['TEC1', None, 'STE12', None])

    gene_name = f_name.split('.')[0]

    plt.savefig(RESULT_DIR / 'promoter_tf_response' / 'figs' / f'{gene_name}.pdf')
    plt.clf()
