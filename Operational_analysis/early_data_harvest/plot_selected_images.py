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

directory='/Users/lindamegner/MATS/retrieval/FlightData/221109_first_images/RacFiles_out_72/'

calibration_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'
channels=['IR1','IR2','IR3','IR4','UV1','UV2']

instrument = Instrument(calibration_file)

calibrate=True  

CCDitems = read_CCDitems(directory)  # read in data
print('Total number of CCDitems: ',len(CCDitems))


if calibrate:
    calibrate_CCDitems(CCDitems, instrument) #calibrate



plotdir=directory[:-16]+'plot_dir'+directory[-4:-1]+'_cal/'

#shutil.rmtree(plotdir)

os.mkdir(plotdir)
os.mkdir(plotdir+'/single_images')

full_data=True #Set to false for movies

#Select and plot data
#key_value_dict = {'channel': 'IR1'}
#key_value_dict = {'channel': 'IR1', 'channel': '',}


climfact= {'IR1': 2.5, 'IR2': 5.,'IR3': 0.5,'IR4': 0.5,'UV1': .15,'UV2': 2.5,'NADIR': 2.,}
# orbit 74 climfact= {'IR1': 2., 'IR2': 4.,'IR3': 0.5,'IR4': 0.5,'IR1': 2,'UV1': .15,'UV2': 2.5,'NADIR': 2.,}
#climfact= {'IR1': .075, 'IR2': .13,'IR3': .02,'IR4': .02,'UV1': .1,'UV2': 0.1,'NADIR': 0.1,}
climfact_cal= {'IR1': 0.8, 'IR2': 0.8,'IR3': 0.6,'IR4': 0.6,'UV1': 0.8,'UV2': 0.8,'NADIR': 1.,}

texpmslist=[1000,2000,3000,4000, 5000, 6000,7000, 8000, 9000, 10000, 12000, 
            13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000 ]

for channel in channels:
    CCDitems_select=select_CCDitems(CCDitems, 'channel',channel)

    dirname=plotdir+'single_images/'+channel
    os.mkdir(dirname)
    for texpms in texpmslist:
        CCDitems_select_texpms=select_CCDitems(CCDitems_select, 'TEXPMS',texpms)
        if len(CCDitems_select_texpms)>0:
            texpmsdir=dirname+'/texpms_'+str(texpms)
            os.mkdir(texpmsdir)
            for CCDitem in CCDitems_select_texpms:
                if full_data: #set to false for movies
                    fig , ax= plt.subplots(squeeze=False, frameon=False)
                    ax = plt.subplot2grid((2, 1), (0, 0))    
                else:    
                    fig , ax= plt.subplots(1, 1, frameon=False, figsize=(4, 1))

                title='Calibrated signal in $10^{10}$ photons nm$^{-1}$ cm$^{-2}$ s$^{-1}$ str$^{-1}$'
                CCDitem['channel']+'_'+str(CCDitem['TEXPMS'])+'_'+CCDitem['WDW InputDataWindow']
                totbin=CCDitem['NCBIN CCDColumns']*CCDitem['NCBIN FPGAColumns']*CCDitem['NRBIN']
                if totbin!=1:
                    if calibrate:
                        clim=[0,climfact_cal[channel]*CCDitem['TEXPMS']]
                    else:
                        clim=[0,(climfact[channel]*CCDitem['TEXPMS']+200)]
                else:
                    clim=999  
                sp=plot_CCDimage(CCDitem['image_calibrated'], fig, ax, title=title, clim=clim)
                        
                if full_data: #set to false for movies
                    info=['id','channel', 'EXP Date','File',
                          'JPEGQ', 'NCBIN CCDColumns', 
                          'NCOL', 'NRBIN' , 'NROW', 'TBLNK', 'TEMP', 
                          'temperature','temperature_ADC','temperature_HTR', 'TEXPMS', 
                          'WDW InputDataWindow','WDW Mode']
                    
                    if CCDitem['NCBIN CCDColumns']>1:
                        xpos=-30
                        xdiff=-20
                    else:
                        xpos=-70
                        xdiff=-40
                    for ikey in info:
                        xpos=xpos+xdiff
                        ax.text(0, xpos, ikey+': '+str(CCDitem[ikey]) )
               
                    #ax.text(0, 20, 'climax: '+str(clim) )
                fig.savefig(texpmsdir+'/image_'+CCDitem['id']+'.png', dpi=600)
                plt.close(fig)
