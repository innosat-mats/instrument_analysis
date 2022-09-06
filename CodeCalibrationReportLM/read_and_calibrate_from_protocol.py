#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""




from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage, read_all_files_in_protocol, plotCCDitem
from mats_l1_processing.read_in_functions import readprotocol
import matplotlib.pyplot as plt
from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.instrument import Instrument





#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_0C/'
#protocol='flatfields_200330_SigMod1_LMprotocol.txt'

directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/'
protocol='protocol_dark_bright_100um_incl_IR3.txt'

#directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/BinningFlatfieldsIR3_210910/'
#protocol='readprotocol_IR3.txt'


read_from="rac" 
df_protocol=readprotocol(directory+protocol)

df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)

calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'
instrument = Instrument(calibration_file)


calibrate=False


for CCDitem in CCDitems[:]:
     
    if calibrate:
 
        image_lsb,image_bias_sub,image_desmeared, image_dark_sub, image_flatf_comp =L1_calibrate(CCDitem, instrument)

        fig,ax=plt.subplots(5,1)
        plot_CCDimage(image_lsb,fig, ax[0], 'Original LSB')    
        plot_CCDimage(image_bias_sub,fig, ax[1], 'Bias subtracted')  
        plot_CCDimage(image_desmeared,fig, ax[2],' Desmeared LSB')  
        plot_CCDimage(image_dark_sub,fig, ax[3], ' Dark current subtracted LSB')  
        plot_CCDimage(image_flatf_comp,fig, ax[4], ' Flat field compensated LSB')         
        fig.suptitle(CCDitem['channel'])

    else:    
    
        fig=plt.figure()
        ax=fig.gca()
        plotCCDitem(CCDitem,fig, ax,CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000)+'s')   
