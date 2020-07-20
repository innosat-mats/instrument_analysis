#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 11:53:02 2020

@author: bjorn
"""

## INVESTIGATE TBLANK; LBLANK VS SIGNAL




import sys
sys.path.append('/Users/bjorn/Documents/PhD/MATS/calibration/MATS-L0-processing')
sys.path.append('/Users/bjorn/Documents/PhD/MATS/calibration/read_images')

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from LindasCalibrationFunctions import plotCCDitem
from LindasCalibrationFunctions import plot_simple
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics
from L1_calibration_functions import get_true_image_remove, desmear_true_image_remove, desmear_true_image
import copy

cal_day = '/080720_binning_nonlinearity'
#cal_day = '/20200713BinningTest'

ref_time = 4
#list_exposure = 'exposure.txt'
list_exposure = 'exposure_11.txt'


####################################################
############# Load Reference Images ################
####################################################

dirname = ('/Users/bjorn/Documents/PhD/MATS/calibration'+ cal_day 
                + '/PayloadImages')

os.chdir(dirname)

list_names = ['column.txt', 'row.txt']
#list_names = ['column_11.txt','row_11.txt', 'fpga.txt']
#list_names = ['column.txt']
plt.figure()
plt.grid()
for list_name in list_names:
    
    CCDitems = []
    binned = []
    
    itemlist = np.genfromtxt(list_name, dtype='str')

          
    for line in itemlist:
        
        CCD = read_CCDitem_from_imgview(dirname, line)
            
        CCDitems.append(copy.deepcopy(CCD))
     
    CCDd_list = copy.deepcopy(CCDitems[0::2])
    CCDs_list = copy.deepcopy(CCDitems[1::2])
    
    
    
    for i in range(0, len(CCDd_list)-2):
        if list_name == 'row.txt':
            plt.scatter(CCDd_list[i]['NRBIN'],CCDd_list[i]['LBLNK'], color='blue', marker='x')
            plt.scatter(CCDs_list[i]['NRBIN'],CCDs_list[i]['LBLNK'], color='blue', marker='o')
            #plt.scatter(CCDs_list[i]['NRBIN'],CCDs_list[i]['TBLNK']-CCDd_list[i]['TBLNK'], color='blue', marker='o')
        else:
            # plt.scatter(CCDd_list[i]['NColBinCCD'],CCDd_list[i]['TBLNK'], color='red', marker='x')
            # plt.scatter(CCDs_list[i]['NColBinCCD'],CCDs_list[i]['TBLNK'], color='red', marker='o')
            
            plt.scatter(CCDd_list[i]['NColBinCCD'],CCDd_list[i]['LBLNK'], color='red', marker='x')
            plt.scatter(CCDs_list[i]['NColBinCCD'],CCDs_list[i]['LBLNK'], color='red', marker='o')
            #plt.scatter(CCDs_list[i]['NColBinCCD'],CCDs_list[i]['TBLNK']-CCDd_list[i]['TBLNK'], color='red', marker='o')
    
    
    # for i in range(0, len(CCDd_list)-2):
    #     if list_name == 'row.txt':
    #         plt.scatter(CCDd_list[i]['IMAGE'].mean(),CCDd_list[i]['LBLNK'], color='blue', marker='x')
    #         plt.scatter(CCDs_list[i]['IMAGE'].mean(),CCDs_list[i]['LBLNK'], color='blue', marker='o')
    #     else:
    #         plt.scatter(CCDd_list[i]['IMAGE'].mean(),CCDd_list[i]['LBLNK'], color='red', marker='x')
    #         plt.scatter(CCDs_list[i]['IMAGE'].mean(),CCDs_list[i]['LBLNK'], color='red', marker='o')
        
    
    # 19 mer varje pixel

#plt.xlim([0,20000])
plt.title('LBLNK vs. binning size ; row bin (blue), col bin (red)')

    # if list_name == 'column.txt':
    #     factors.append(CCDs_list[i]['LBLNK'])
        
    # if list_name == 'row.txt':


   