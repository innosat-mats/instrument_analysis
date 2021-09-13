#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""


import sys


from mats_l1_processing.LindasCalibrationFunctions import (
    plotCCDitem,
    read_all_files_in_protocol,
)


from mats_l1_processing.read_in_functions import add_temperature_info_to_CCDitems

import matplotlib.pyplot as plt


directory = "/Users/bjorn/Documents/PhD/MATS/calibration/DarkMeas_20200424/"

protocol = "protocol.txt"


read_from = "imgview"

df, CCDitems = read_all_files_in_protocol(protocol, read_from, directory)


CCDitems = add_temperature_info_to_CCDitems(CCDitems, read_from, directory)


fig, ax = plt.subplots(1, 1)

plotCCDitem(CCDitems[0], fig, ax)
