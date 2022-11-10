#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: lindamegner


This is more or less a duplicate of read_and_calibrate_all_files_in_directory but it is meant to be used as a script.
"""
#%%
from plotting.sort_images import sort_images_in_dirs, sort_images_plot
from mats_l1_processing.read_in_functions import read_CCDitems
import matplotlib.pyplot as plt

from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage
from mats_l1_processing.instrument import Instrument
import numpy as np

def bin_image(image,nrbin,ncolbin):

    nrow,ncol = image.shape()
    nrow_binned = np.floor(nrow/nrbin)
    ncol_binned = np.floor(ncol/ncolbin)

    colbin = np.zeros([nrow,ncol_binned])
          
    for j in range(0,ncol_binned):
        colbin[:,j] = image[:,j*ncolbin:(j+1)*ncolbin].sum(axis=1)

    # declare zero array for row binning 
    binned = np.zeros([nrow_binned,ncol_binned])
    
    for j in range(0,nrow_binned):
        binned[j,:] = colbin[j*nrbin:(j+1)*nrbin,:].sum(axis=0)

    return binned


def calibrate_CCDitems(CCDitems,instrument, plot=False):
    """
    Calibrate all CCDitems in the list

    Parameters
    ----------
    CCDitems : List of dictionaries
        Contains images and housing data
    instrument: instrument object, see mats_l1_processing.instrument
    plot : logical, optional
        If true the calibrations steps are plotted. The default is False.

    Returns
    -------
    Does not return anything but CCDitems will now contain calibrated images
    """
    for CCDitem in CCDitems:
        (
            image_lsb,
            image_bias_sub,
            image_desmeared,
            image_dark_sub,
            image_calib_nonflipped,
            image_calibrated,
            errors
        ) = L1_calibrate(CCDitem, instrument)

        if plot:
            fig, ax = plt.subplots(5, 1)
            plot_CCDimage(image_lsb, fig, ax[0], "Original LSB")
            plot_CCDimage(image_bias_sub, fig, ax[1], "Bias subtracted")
            plot_CCDimage(image_desmeared, fig, ax[2], " Desmeared LSB")
            plot_CCDimage(
                image_dark_sub, fig, ax[3], " Dark current subtracted LSB"
            )
            plot_CCDimage(
                image_calib_nonflipped, fig, ax[4], " Flat field compensated LSB"
            )
            fig.suptitle(CCDitem["channel"])
    



directory='out/'

calibration_file='calibration_data_linda.toml'

instrument = Instrument(calibration_file)

calibrate=True  

CCDitems = read_CCDitems(directory)  # read in data
print('Total number of CCDitems: ',len(CCDitems))

if calibrate:
    calibrate_CCDitems(CCDitems, instrument) #calibrate

binned = bin_image(CCDitems[0]["image_calibrated"],2,2)

# %%
fig, ax = plt.subplots(1, 1,frameon=False)
ax.set_aspect('equal')
ax.axis('off')
ax.pcolor(CCDitems[0]["image_calibrated"],cmap='gray',clim=[-1,10])
fig.savefig('IR1.png', format='png', dpi=1200)
# %%
