import pandas as pd
from matplotlib import pyplot as plt
from promo_seq_anal.constants import *

data = pd.read_csv(RESULT_DIR / 'promoter_tf_response' / 'TEC1.tsv', delimiter='\t', index_col=False)
data_ste = pd.read_csv(RESULT_DIR / 'promoter_tf_response' / 'STE12_consensus.tsv', delimiter='\t', index_col=False)

plt.figure()
plt.hist(data['max_bind'], 100)
plt.title('TEC1 max-binding histogram')
plt.xlabel('Maximum Relative binding affinity')
plt.ylabel('# Promoters')

plt.figure()
plt.hist(data['sum_bind'], 100)
plt.title('TEC1 sum-binding histogram')
plt.xlabel('Relative binding affinity sum')
plt.ylabel('# Promoters')


plt.figure()
plt.hist(data_ste['max_bind'], 100)
plt.title('STE12 max-binding histogram')
plt.xlabel('Maximum Relative binding affinity')
plt.ylabel('# Promoters')

plt.figure()
plt.hist(data_ste['sum_bind'], 100)
plt.title('STE12 sum-binding histogram')
plt.xlabel('Relative binding affinity sum')
plt.ylabel('# Promoters')


plt.show()

print()
