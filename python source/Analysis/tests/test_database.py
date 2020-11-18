import unittest
from functionality.database import *
import numpy as np


class TestSubscriptions(unittest.TestCase):
    """Tests all functions in the module ´database´"""

    def test_load_data(self):
        """Tests the ´load_data´ function"""
        example_file = "2010.alpha_time.pcl"
        load_data(example_file)

    def test_save_list_of_names(self):
        """Tests the ´load_data´ function"""
        example_file = "2010.alpha_time.pcl"
        save_list_of_names(example_file)

    def test_show_decreasing_yorfs_time(self):
        """Tests the ´show_decreasing_yorfs_time´ function"""
        example_file = "2010.alpha_time.pcl"
        data = load_data(example_file)
        show_decreasing_yorfs_time(data)

    def test_show_decreasing_yorfs_concentrations(self):
        """Tests the ´show_decreasing_yorfs_concentrations´ function"""
        example_file = "2010.alpha_conc.pcl"
        data = load_data(example_file)
        show_decreasing_yorfs_concentrations(data)

    def test_cpu_test(self):
        cpu_test()
