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
from mats_l1_processing.L1b_calibration_functions import shift_image


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
fig3, ax3 = plt.subplots(6, 1)


bottom_side_covered_start=13
bottom_side_covered_stop=18+1
left_side_covered_start=19
left_side_covered_stop=22+1
deuterium_LFOV_start=23
deuterium_LFOV_stop=28+1
deuterium_RFOV_start=29
deuterium_RFOV_stop=34+1
laser1_start=35
laser1_stop=40+1
laser2_start=41
laser2_stop=48+1




start=deuterium_LFOV_start
stop=deuterium_LFOV_stop

# ind=-1
# for i in [23, 28]:
#     CCDitem=CCDitems[i]
#     ind=ind+1

for ind, CCDitem in enumerate(CCDitems[left_side_covered_start:left_side_covered_stop]):
    (
        image_lsb,
        image_bias_sub,
        image_desmeared,
        image_dark_sub,
        image_calib_nonflipped,
        image_calibrated,
        errors
    ) = L1_calibrate(CCDitem, instrument)
    #Shift image, i.e. put image on common field of view
    image_common_fov, error_flags_flipnshift = shift_image(CCDitem, image_calibrated)
    

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
    plot_CCDimage(image_common_fov, fig3, ax3[ind], title=str(start)+CCDitem['channel'], clim=fix_clim)  
    fig3.savefig('images/Alignment_all_'+str(start)+'_'+CCDitem['channel']+'.jpg')   
          
    sp=plot_CCDimage_transp(image_common_fov, fig2, ax2, title='All channels'+CCDitem['channel'], clim=fix_clim, alpha=0.3)    
    plt.show()
    fig2.savefig('images/Alignment_together_'+str(start)+'_'+CCDitem['channel']+'.jpg') 
    

#fig.savefig('Testfile.jpg')    
