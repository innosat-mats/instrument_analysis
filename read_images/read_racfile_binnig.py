"""
Created on Fri Aug 14 15:15:36 2020

@author: olemar
"""
# fmt: off
import sys

sys.path.insert(1, "/home/olemar/Projects/MATS/MATS-L1-processsing")
import read_in_functions

# fmt: on


CCD_image_data = read_in_functions.read_CCDitems(
    "/home/olemar/Projects/MATS/MATS-data/binning_test_20200812_racfiles/binning/"
)
