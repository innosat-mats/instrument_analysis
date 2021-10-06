#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""



from LindasCalibrationFunctions import  plotCCDitem



from read_in_functions import read_all_files_in_directory
import matplotlib.pyplot as plt
from L1_calibrate import L1_calibrate
from LindasCalibrationFunctions import plot_CCDimage



directory='/Users/lindamegner/MATS/retrieval/Calibration/CoolingTests/Rac22C/'

directory='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_AIT/20190520/'
directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210616_Baffletests_and_LimbHouseLightLeakage/'
#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/'
read_from='rac'  
CCDitems=read_all_files_in_directory(read_from,directory)



calibrate=True
plot=True

if calibrate:
    for CCDitem in CCDitems[:3]:
        image_lsb,image_bias_sub,image_desmeared, image_dark_sub, image_flatf_comp =L1_calibrate(CCDitem)

        if plot==True:
            fig,ax=plt.subplots(5,1)
            plot_CCDimage(image_lsb,fig, ax[0], 'Original LSB')    
            plot_CCDimage(image_bias_sub,fig, ax[1], 'Bias subtracted')  
            plot_CCDimage(image_desmeared,fig, ax[2],' Desmeared LSB')  
            plot_CCDimage(image_dark_sub,fig, ax[3], ' Dark current subtracted LSB')  
            plot_CCDimage(image_flatf_comp,fig, ax[4], ' Flat field compensated LSB')         
            fig.suptitle(CCDitem['channel'])

else:    

    for CCDitem in CCDitems[:]:
        fig=plt.figure()
        ax=fig.gca()
        plotCCDitem(CCDitem,fig, ax, title=CCDitem['channel'])
 