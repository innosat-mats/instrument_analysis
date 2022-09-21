#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:16:26 2020

@author: bjorn
"""


import sys
sys.path.append('/Users/bjorn/Documents/PhD/MATS/calibration/MATS-L0-processing')

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from experimental_utils import  plotCCDitem
from experimental_utils import  plot_simple
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics
from L1_calibration_functions import get_true_image_remove, desmear_true_image_remove

####################################################
############# Load Reference Images ################
####################################################

# Location of measurements
#cal_day = '/19052020_binning/'
#cal_day = '/27052020_binning/'
cal_day = '/03062020_binning/'
meas = 'binning_simple'

# number of tests
num_test = 5

timing1 = np.zeros([5,6])
timing2 = np.zeros([5,6])
timing3 = np.zeros([5,6])

# declare empty lists
IDstrings = []
CCDitems= []

for test in range(0,5):
    
    
    dirname = ('/Users/bjorn/Documents/PhD/MATS/calibration'+ cal_day 
                + meas + '/PayloadImages/'+'test_' + str(test+1))
    
    os.chdir(dirname)
    
    for file in glob.glob("*.pnm"):
        IDstrings.append(file.strip('.pnm'))
        
    IDstrings = sorted(IDstrings, key=lambda x:x[-1:])
    
    for IDstring in IDstrings:
        CCD = read_CCDitem_from_imgview(dirname, IDstring)
        
        CCDitems.append(CCD)

    ####################################################
    ############ Iterate through Channels ##############
    ####################################################
    
    
    for i in range(0,6):
        
        timing1[test,i] = CCDitems[i]['TIMING1']
        timing2[test,i] = CCDitems[i]['TIMING2']
        timing3[test,i] = CCDitems[i]['TIMING3']
    
    IDstrings = []
    CCDitems= []
    
# TIMING 1; 
# Time before R goes high (Time0) / Stop clamping the signal (ClampStop)
# Time before R goes low (Time1) / Start clamping the signal (ClampStart)

c, f= divmod(timing1, 256)

# Time0 = 19; Time1 = 3 
    
# TIMING 2; 
# Time before R3 goes high (Time2) / Start signal conversion (ConvertStrobe)
# Time before R1 goes high (Time3) / Setting to replace Time0, Time2, Time3 and 
#                            Time4 for pixels with fast-clocking (TimeFast)

c, f= divmod(timing2, 256)

# Time2 = 9,; Time3 = 3
# ConvertStrobe = 6; TimeFast  = 2

# TIMING 3; 
# Time before R2 goes high (Time4) / Setting to replace Time1 for pixels with 
#                                       fast-clocking (TimeFastR)
# Overlap time of R clocks (TimeOverlap) / TimeofI clocks during row shifting 
#                                           (RowShift)

c, f= divmod(timing3, 256)

# Time4 = 3; TimeOverlap = 1
# TimeFastR = 3; RowShift = 164

# The alternating data does work. The values are however exactly the same as 
# the defaults.