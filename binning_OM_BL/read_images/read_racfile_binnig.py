"""
Created on Fri Aug 14 15:15:36 2020

@author: olemar
"""
# fmt: off
import sys

from mats_l1_processing import read_in_functions

# fmt: on


CCD_image_data = read_in_functions.read_CCDitems(
    "/home/olemar/Projects/MATS/MATS-data/binning_test_20200812_racfiles/binning/"
)
