from promo_seq_anal.pfm_reader import read_pfm
from pathlib import Path
import os

from matplotlib import pyplot as plt


binding_seq = read_pfm('data/STE12.pfm', ignore_sum=True)

print()
