#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: lindamegner


This is more or less a duplicate of read_and_calibrate_all_files_in_directory but it is meant to be used as a script.
"""

from mats_l1_processing.experimental_utils import plotCCDitem
from sort_images import sort_images_in_dirs, sort_images_plot
from mats_l1_processing.read_in_functions import read_CCDitems
import matplotlib.pyplot as plt

from mats_l1_processing.L1_calibrate import L1_calibrate
from mats_l1_processing.experimental_utils import plot_CCDimage, calibrate_CCDitems
from mats_l1_processing.instrument import Instrument
import pickle
import os
import shutil



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
    



#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210215OHBLimbImage/RacFiles_out/'
#directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/RacFiles_210906-210910/RacFiles_out/'

#aligning images
#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/210226OHB_IR3Image/RacFiles_out/'

#binned images with nskip
#directory='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/testdata/binning_test_20200812/RacFiles_out/'
#directory='/Users/lindamegner/MATS/retrieval/Calibration/FinalFinalSept2021/LaserAlignmentTest_RacFiles0909to0910Duplicated/RacFiles_out/'

directory='/Users/lindamegner/MATS/retrieval/FlightData/221109_first_images/RacFiles_out/'

calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'
channels=['IR1', 'IR2','IR3', 'IR4','UV1','UV2','NADIR']

instrument = Instrument(calibration_file)

calibrate=True  

CCDitems = read_CCDitems(directory)  # read in data
print('Total number of CCDitems: ',len(CCDitems))

if calibrate:
    calibrate_CCDitems(CCDitems, instrument) #calibrate





#Select and plot data
key_value_dict = {'channel': 'IR1'}
#key_value_dict = {'channel': 'IR4', 'TEXPMS': 2000}



plotdir=directory[:-13]+'plot_dir'

shutil.rmtree(plotdir)

os.mkdir(plotdir)
os.mkdir(plotdir+'/single_images')
#sort_images_in_dirs(CCDitems, key_value_dict, path=plotdir+'/single_images', whattoplot="image_calibrated")

os.mkdir(plotdir+'/multi_images')

sort_images_plot(CCDitems, key_value_dict,clim=[-1,10],path=plotdir+'/multi_images', whattoplot="image_calibrated")

