import pandas as pd
from .constants import *
from matplotlib import pyplot as plt
import numpy as np
import math


def load_data(pcl_file):
    return pd.read_csv(DATA_ROOT + pcl_file, delim_whitespace=True)


def show_decreasing_yorfs_time(data):
    # TODO-refactor: create list with time keys
    decreasing_data = data[(data["0-min"] > data["15-min"]) &
                           (data["15-min"] > data["30-min"]) &
                           (data["30-min"] > data["45-min"]) &
                           (data["45-min"] > data["60-min"]) &
                           (data["60-min"] > data["90-min"]) &
                           (data["90-min"] > data["120-min"])]

    legend = []
    yorfs = []
    time_vector = [0, 15, 30, 45, 60, 90, 120]

    for index, row in decreasing_data.iterrows():
        if type(row["NAME"]) == str:
            plt.plot(time_vector,
                     [row["0-min"], row["15-min"], row["30-min"], row["45-min"], row["60-min"], row["90-min"],
                      row["120-min"]])
            legend.append(row["NAME"])
        yorfs.append(row["YORF"])

    plt.ylabel("Log2(Expression Ratio)")
    plt.xlabel("Time (Minutes)")
    plt.legend(legend)
    plt.title("Alpha-factor response: Monotonously decreasing expression ORFs")

    return legend, yorfs


def show_decreasing_yorfs_concentrations(data):
    # TODO-refactor: create list with time keys
    decreasing_data = data[(data["0.15-nM-aF"] > data["0.5-nM-aF"]) &
                           (data["0.5-nM-aF"] > data["1.5-nM"]) &
                           (data["1.5-nM"] > data["5-nM-aF"]) &
                           (data["5-nM-aF"] > data["15.8-nM-aF"]) &
                           (data["15.8-nM-aF"] > data["50-nM-aF"]) &
                           (data["50-nM-aF"] > data["158-nM-aF"]) &
                           (data["158-nM-aF"] > data["500-nM-aF"])]

    legend = []
    yorfs = []
    concetrations_vector = [15, 50, 150, 500, 15800, 50000, 158000, 500000]

    for index, row in decreasing_data.iterrows():
        if type(row["NAME"]) == str:
            plt.plot(concetrations_vector,
                     [row["0.15-nM-aF"], row["0.5-nM-aF"], row["1.5-nM"], row["5-nM-aF"], row["15.8-nM-aF"],
                      row["50-nM-aF"], row["158-nM-aF"], row["500-nM-aF"]])
            legend.append(row["NAME"])
        yorfs.append(row["YORF"])

    plt.ylabel("Log2(Expression Ratio)")
    plt.xlabel("Concentration (pM)")
    plt.legend(legend)
    plt.title("Alpha-factor response: Monotonously decreasing expression ORFs")
    plt.xscale('log')

    return legend, yorfs


def save_list_of_names(pcl_file):
    data = load_data(pcl_file)
    yorf_names = " ".join(data["YORF"])
    with open(DATA_ROOT + "yorf_names.txt", 'w') as out_file:
        out_file.write(yorf_names)


def cpu_test():
    n = 10000
    m = 10000

    a = np.random.rand(m, m)
    b = np.random.rand(m, m)

    for i in range(n):
        c = np.dot(a, b)

    print("fin")
