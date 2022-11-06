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
from PIL import Image
import warnings


def plot_CCDimage_transp(image, fig, axis, title="", clim=999, aspect="auto", alpha=1):
    sp = axis.imshow(image, cmap="viridis", origin="lower", interpolation="none", alpha=alpha)
    # sp=axis.pcolormesh(image, , cmap='viridis')
    if clim == 999:
        mean = image.mean()
        std = image.std()
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    #fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    axis.set_aspect(aspect)
    return sp



#Setup directories 
#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/'
#df_protocol=readprotocol(directory+'protocol_dark_bright_incl_IR3.txt')
directory='/Users/lindamegner/MATS/retrieval/Calibration/2022_Tests_At_Launch_Site/LimbTests221012/'
df_protocol=readprotocol(directory+'protocol.txt')
read_from='rac'  
df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)


if directory=='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/':
    warnings.warn("Warning: IR3 images were taken seperately, to align them use x_pos=37+25,y_pos=66+60 in flip_and_shift")
#Select channel and choos picture for that channel
#channel='UV2'

#CCDitem_select= list(filter(lambda x: ( x['channel']==channel),CCDitems))
#CCDitem=CCDitem_select[0]

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
    
    
    

#fig.savefig('Testfile.jpg')    
