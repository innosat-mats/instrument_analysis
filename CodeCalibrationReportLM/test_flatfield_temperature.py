  #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 21:00:16 2022

@author: lindamegner


"""


from mats_l1_processing.experimental_utils import plot_CCDimage,read_all_files_in_protocol
from  mats_l1_processing.experimental_utils import readprotocol
from mats_l1_processing.instrument import Instrument
from mats_l1_processing.L1_calibrate import L1_calibrate
import matplotlib.pyplot as plt
import numpy as np
#from jpeglib import read12bit_jpegfile

from mats_l1_processing.items_units_functions import read_files_in_protocol_as_ItemsUnits



directory_22C='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200506_flatfields_roomtemp/'
directory_8C='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200429_flatfields_8C/'
directory_0C='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_0C/'

directories=[directory_0C, directory_8C, directory_22C]
temperature=[0, 8, 22]

read_from='rac'  
    

exptime=6000
l_ir1=[]
l_ir2=[]
l_ir3=[]
l_ir4=[]
l_uv1=[]
l_uv2=[]
for directory in directories[1:]:
    df_protocol=readprotocol(directory+'protocol.txt')
    CCDunits=read_files_in_protocol_as_ItemsUnits(df_protocol, directory, 3, read_from)
    for CCDunit in CCDunits:
        if CCDunit.imageItem['TEXPMS']==exptime:
            if CCDunit.imageItem['channel']=='IR1':
                l_ir1.append(CCDunit.subpic.mean())
            elif CCDunit.imageItem['channel']=='IR2':
                l_ir2.append(CCDunit.subpic.mean())
            elif CCDunit.imageItem['channel']=='IR3':
                l_ir3.append(CCDunit.subpic.mean())
            elif CCDunit.imageItem['channel']=='IR4':
                l_ir4.append(CCDunit.subpic.mean())
            elif CCDunit.imageItem['channel']=='UV1':
                l_uv1.append(CCDunit.subpic.mean())
            elif CCDunit.imageItem['channel']=='UV2':
                l_uv2.append(CCDunit.subpic.mean())                
    
a_ir1 = np.asarray(l_ir1)
a_ir2 = np.asarray(l_ir2)
a_ir3 = np.asarray(l_ir3)
a_ir4 = np.asarray(l_ir4)
a_uv1 = np.asarray(l_uv1)
a_uv2 = np.asarray(l_uv2)
    

fig1, ax1 = plt.subplots(1, 2)
ax1[0].plot(temperature, a_ir1, '*')
ax1[1].plot(temperature, a_ir2, '*')
# ax1.plot(temperature, a_ir3, '*')
# ax1.plot(temperature, a_ir4, '*')
# ax1.plot(temperature, a_uv1, '*')
# ax1.plot(temperature, a_uv2, '*')

"""
    df_bright=df_protocol[df_protocol.DarkBright=='B']
    CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)


calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'


plot=True

#CCDitems = read_CCDitems(directory)  # read in data

instrument = Instrument(calibration_file)



fig1, ax1 = plt.subplots(6, 1)
fig2, ax2 = plt.subplots(1)

for ind, CCDitem in enumerate(CCDitems[:6]):
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated,
        image_common_fov,
        errors
    ) = L1_calibrate(CCDitem, instrument)

    if plot:
        fig, ax = plt.subplots(6, 1)
        plot_CCDimage(image_lsb, fig, ax[0], "Original LSB")
        plot_CCDimage(image_bias_sub, fig, ax[1], "Bias subtracted")
        plot_CCDimage(image_desmeared, fig, ax[2], " Desmeared LSB")
        sp=plot_CCDimage(image_dark_sub, fig, ax[3], " Dark current subtracted LSB")
        fix_clim=sp.get_clim()
        plot_CCDimage(image_calib_nonflipped, fig, ax[4], " Flat field compensated LSB", clim=fix_clim)
        plot_CCDimage(image_common_fov, fig, ax[5], " Image common field of view", clim=fix_clim)
        fig.suptitle(CCDitem["channel"])
            
    plot_CCDimage(image_lsb, fig1, ax1[ind], CCDitem['channel'])    
    fig.savefig('images/Alignment_'+CCDitem['channel']+'.jpg')             
    sp=plot_CCDimage_transp(image_common_fov, fig2, ax2, 'All channels', clim=fix_clim, alpha=0.3)    
    
    
    

"""