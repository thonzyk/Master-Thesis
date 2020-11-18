import numpy as np
from functionality.database import *
from matplotlib import pyplot as plt


def debug_1():
    time_file = "2010.alpha_time.pcl"
    time_data = load_data(time_file)
    plt.figure()
    t_labels, t_yorfs = show_decreasing_yorfs_time(time_data)
    plt.figure()
    concentrations_file = "2010.alpha_conc.pcl"
    concentrations_data = load_data(concentrations_file)
    c_labels, c_yorfs = show_decreasing_yorfs_concentrations(concentrations_data)

    print(len(t_yorfs) + len(c_yorfs))
    print(len(set(t_yorfs + c_yorfs)))

    print(len(t_yorfs))
    print(len(c_yorfs))

    plt.show()


if __name__ == '__main__':
    debug_1()
