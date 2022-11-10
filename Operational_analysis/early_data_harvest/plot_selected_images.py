#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 13:53:36 2022

@author: lindamegner


This is more or less a duplicate of read_and_calibrate_all_files_in_directory but it is meant to be used as a script.
"""

from mats_l1_processing.experimental_utils import plotCCDitem
from sort_images import sort_images_in_dirs, sort_images_plot, select_CCDitems, select_CCDitems_using_list, plot_CCDitems
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

directory='/Users/lindamegner/MATS/retrieval/FlightData/221109_first_images/RacFiles_out_75/'

calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'
channels=['IR1','IR2','IR3','IR4','UV1','UV2','NADIR']

instrument = Instrument(calibration_file)

calibrate=False  

CCDitems = read_CCDitems(directory)  # read in data
print('Total number of CCDitems: ',len(CCDitems))


if calibrate:
    calibrate_CCDitems(CCDitems, instrument) #calibrate



plotdir=directory[:-16]+'plot_dir'+directory[-4:]

shutil.rmtree(plotdir)

os.mkdir(plotdir)
os.mkdir(plotdir+'/single_images')



#Select and plot data
#key_value_dict = {'channel': 'IR1'}
#key_value_dict = {'channel': 'IR1', 'channel': '',}




climfact= {'IR1': 2., 'IR2': 4.,'IR3': 0.5,'IR4': 0.5,'IR1': 2,'UV1': 2.,'NADIR': 2.,}


for channel in channels:
    CCDitems_select=select_CCDitems(CCDitems, 'channel',channel)
    #CCDitems_select=select_CCDitems(CCDitems_select, 'WDW InputDataWindow','15..0')
    



          
#texpmslist=[1000, 3000, 6000, 9000, 12000]
#winmode=[4, 7, 128]
#select_CCDitems_using_list(CCDitems, 'TEXPMS', texpmslist)
    #CCDitems_select=select_CCDitems(CCDitems,key,value)
    dirname=plotdir+'single_images/'+channel
    os.mkdir(dirname)
    for CCDitem in CCDitems_select:
        #fig , ax= plt.subplots(1, 1, squeeze=False, frameon=False)
        fig , ax= plt.subplots(1, 1, frameon=False, figsize=(4, 1))
        
        
        #ax = plt.subplot2grid((2, 1), (0, 0))
        
        
        #fig=plt.figure(figsize=(4, 1))
        #ax = fig.gca()
        title=CCDitem['channel']+'_'+str(CCDitem['TEXPMS'])+'_'+CCDitem['WDW InputDataWindow']
        clim=[0,(climfact[channel]*CCDitem['TEXPMS']+200)]
        sp=plot_CCDimage(CCDitem['IMAGE'], fig, ax, title=title, clim=clim)
        #fig.suptitle(channel)
        #sp.set(visible='False')        
        #attempt to set aspect ratio to 1
        #ax.set_aspect(1)
        
        #fig.delaxes(ax.flatten()[1])

        plt.close(fig)
        
        
        
        
        # info=['id','channel', 'EXP Date','File',
        #       'JPEGQ', 'NCBIN CCDColumns', 
        #       'NCOL', 'NRBIN' , 'NROW', 'TBLNK', 'TEMP', 
        #       'temperature','temperature_ADC','temperature_HTR', 'TEXPMS', 
        #       'WDW InputDataWindow','WDW Mode']
        
        # if CCDitem['NCBIN CCDColumns']>1:
        #     xpos=-30
        #     xdiff=-20
        # else:
        #     xpos=-70
        #     xdiff=-40
        # for ikey in info:
        #     xpos=xpos+xdiff
        #     ax.text(0, xpos, ikey+': '+str(CCDitem[ikey]) )
            
        fig.savefig(dirname+'/image_'+CCDitem['id']+'.png', dpi=300)
"""


maketree=True
if maketree:
    os.mkdir(plotdir+'/single_images')
    sort_images_in_dirs(CCDitems, key_value_dict, path=plotdir+'/single_images', whattoplot="IMAGE")

makeplot=False
if makeplot:
    os.mkdir(plotdir+'/multi_images')
    sort_images_plot(CCDitems, key_value_dict,path=plotdir+'/multi_images', whattoplot="IMAGE")

"""