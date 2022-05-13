from promo_seq_anal.pfm_reader import read_pfm
from pathlib import Path
import os

from matplotlib import pyplot as plt

SRC_DIR = Path(r'C:\Users\tomas\My Drive\Master thesis\Data\Yeast TF PFMs\Tec1')

files = os.listdir(SRC_DIR)

for f in files:
    pfm = read_pfm(SRC_DIR / f)

    plt.figure()
    plt.plot(pfm.transpose(), alpha=0.5)
    plt.legend(['A', 'T', 'G', 'C'])

plt.show()
