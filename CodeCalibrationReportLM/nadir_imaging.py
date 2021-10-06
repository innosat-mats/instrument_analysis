#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""



#from LindasCalibrationFunctions import  plotCCDitem


from  mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage,plotCCDitem, read_all_files_in_protocol 
from  mats_l1_processing.read_in_functions import readprotocol
import matplotlib.pyplot as plt
from mats_l1_processing.L1_calibrate import L1_calibrate



# def plotCCDitem(CCDitem, fig, axis, title="", clim=999):

#     pic = CCDitem["IMAGE"]
#     sp = axis.pcolormesh(pic, cmap=plt.cm.jet)
#     axis.set_title(title)
#     if clim == 999:
#         mean = pic.mean()
#         std = pic.std()
#         sp.set_clim([mean - 2 * std, mean + 2 * std])
#     else:
#         sp.set_clim(clim)

#     fig.colorbar(sp, ax=axis)
    
#     return sp

# def plotCCDitem(CCDitem, fig, axis, title="", clim=999):

#     import matplotlib.pyplot as plt

#     pic = CCDitem["IMAGE"]

#     sp= plt.imshow(pic, cmap=plt.cm.gray)
#     if clim == 999:
#         mean = pic.mean()
#         std = pic.std()
#         plt.clim([mean - 4 * std, mean + 4 * std])
#     else:
#         plt.clim(clim)
#     #plt.colorbar()



directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210216OHBNadirImage/'
protocol='protocol.txt'

read_from='rac'  
df_protocol=readprotocol(directory+protocol)

df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)


calibrate=False

for CCDitem in CCDitems:
    
    if calibrate:
 
        image_lsb,image_bias_sub,image_desmeared, image_dark_sub, image_flatf_comp =L1_calibrate(CCDitem)

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
        plotCCDitem(CCDitem,fig, ax,CCDitem['channel']+' '+str(CCDitem['TEXPMS']/1000)+'s', aspect='equal')   
