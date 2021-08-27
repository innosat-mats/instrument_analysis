#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 16:48:54 2020

@author: bjorn
"""


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
from LindasCalibrationFunctions import  plotCCDitem
from LindasCalibrationFunctions import  plot_simple
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

LBLNKS = np.zeros([5+10,6])
TBLNKS = np.zeros([5+10,6])

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
        
        LBLNKS[test,i] = CCDitems[i]['LBLNK']
        TBLNKS[test,i] = CCDitems[i]['TBLNK']

    
    IDstrings = []
    CCDitems= []
    
    
cal_day = '/27052020_binning/'

    
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
        
        LBLNKS[5+test,i] = CCDitems[i]['LBLNK']
        TBLNKS[5+test,i] = CCDitems[i]['TBLNK']

    
    IDstrings = []
    CCDitems= []
    
    
cal_day = '/19052020_binning/'

    
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
        
        LBLNKS[10+test,i] = CCDitems[i]['LBLNK']
        TBLNKS[10+test,i] = CCDitems[i]['TBLNK']

    
    IDstrings = []
    CCDitems= []