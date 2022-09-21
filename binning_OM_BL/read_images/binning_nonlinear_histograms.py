#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:23:19 2020

@author: bjorn
"""



import sys
sys.path.append('/Users/bjorn/Documents/PhD/MATS/calibration/MATS-L0-processing')
sys.path.append('/Users/bjorn/Documents/PhD/MATS/calibration/read_images')

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from experimental_utils import plotCCDitem
from experimental_utils import plot_simple
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics
from L1_calibration_functions import get_true_image_remove, desmear_true_image_remove, desmear_true_image
import copy


def bin_ref(ref, ccd):
    
        # simple code for binning
        
        nrow, ncol, nrskip, ncskip, nrbin, ncbin = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NCBIN CCDColumns'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NCBIN CCDColumns'])
    
        # reference image that will be binned according to 'ccd' settings
        imgref = ref['IMAGE']
        
        # in case reference image is already a binned image
        ncbin, nrbin = int(ncbin/ncbinr), int(nrbin/nrbinr)

        # images must cover the same ccd section
        if ncskip == ncskipr and nrskip == nrskipr:
        
            colbin = np.zeros([nrowr,ncol])
          
            for j in range(0,ncol):
                colbin[:,j] = imgref[:,j*ncbin:(j+1)*ncbin].sum(axis=1)
    
            # declare zero array for row binning 
            binned = np.zeros([nrow,ncol])
            
            for j in range(0,nrow):
                binned[j,:] = colbin[j*nrbin:(j+1)*nrbin,:].sum(axis=0)
           
            return binned

        else:
            
            sys.exit('Error: images not from the same CCD region.')




def bin_ref_FPGA(ref, ccd):
    
        # simple code for binning 
        
        nrow, ncol, nrskip, ncskip, nrbin, ncbin = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NCBIN FPGAColumns'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NCBIN FPGAColumns'])
    
        # reference image that will be binned according to 'ccd' settings
        imgref = ref['IMAGE']
        
        # in case reference image is already a binned image
        ncbin, nrbin = int(ncbin/ncbinr), int(nrbin/nrbinr)

        # images must cover the same ccd section
        if ncskip == ncskipr and nrskip == nrskipr:
        
            colbin = np.zeros([nrowr,ncol])
          
            for j in range(0,ncol):
                colbin[:,j] = imgref[:,j*ncbin:(j+1)*ncbin].sum(axis=1)
    
            # declare zero array for row binning 
            binned = np.zeros([nrow,ncol])
            
            for j in range(0,nrow):
                binned[j,:] = colbin[j*nrbin:(j+1)*nrbin,:].sum(axis=0)
           
            return binned

        else:
            
            sys.exit('Error: images not from the same CCD region.')
            
            
def img_diff(image1, image2):
    
    return image1-image2



####################################################
#############        LOAD DATA      ################
####################################################


cal_day = '/080720_binning_nonlinearity'
#cal_day = '/20200713BinningTest'

dirname = ('/Users/bjorn/Documents/PhD/MATS/calibration'+ cal_day 
                + '/PayloadImages')

os.chdir(dirname)

CCDitems = []
binned = []


# list containing file names of binning tests
list_name='row.txt'
# list_name='column.txt'
# list_name='fpga.txt'
# list_name='exposure.txt'
# list_name='column_11.txt'

itemlist = np.genfromtxt(list_name, dtype='str')

       
for line in itemlist:
    CCD = read_CCDitem_from_imgview(dirname, line)
        
    CCDitems.append(CCD)

CCDd_list = np.copy(CCDitems[0::2])
CCDs_list = np.copy(CCDitems[1::2])

# plot the raw images
fig, axs = plt.subplots(2,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = 1, wspace=.001)
axs = axs.ravel()
fig.suptitle('images')
  
for i in range(0,len(CCDd_list)):
    im = axs[i].imshow(CCDd_list[i]['IMAGE'],cmap='jet')
    mean=CCDd_list[i]['IMAGE'].mean()
    std=CCDd_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i].set_title('dark: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns']))
    
    im = axs[i+len(CCDd_list)].imshow(CCDs_list[i]['IMAGE'],cmap='jet')
    mean=CCDs_list[i]['IMAGE'].mean()
    std=CCDs_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i+len(CCDd_list)])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i+len(CCDd_list)].set_title('shutter: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns']))

# SUBTRACT DARK IMAGES
CCDs_sub_img = []

for i in range(0,len(CCDd_list)):

    CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))

# plot the images with removed dark current    
fig1, axs1 = plt.subplots(1,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig1.subplots_adjust(hspace = 1, wspace=.001)
axs1 = axs1.ravel()
fig1.suptitle('images with subtracted dark images')

for i in range(0,len(CCDd_list)):
    im = axs1[i].imshow(CCDs_sub_img[i],cmap='jet')
    mean=CCDs_sub_img[i].mean()
    std=CCDs_sub_img[i].std()
    fig1.colorbar(im, ax=axs1[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs1[i].set_title('shutter - dark: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns']))
    


# BINNING

# copy settings from bright image
bin_input = copy.deepcopy(CCDs_list)

# first CCD item is the reference that will be binned manually
ref = copy.deepcopy(CCDs_list[0])
ref['IMAGE'] = CCDs_sub_img[0].copy()

# replace images with the images with subtrated dark
for i in range(0,len(CCDd_list)): 
    bin_input[i]['IMAGE'] = CCDs_sub_img[i].copy()
    

# create manually binned images
for i in range(0,len(CCDd_list)): 
    
    if list_name == 'fpga.txt':
        binned.append(bin_ref_FPGA(copy.deepcopy(ref), bin_input[i].copy()))
    
    else:
        binned.append(bin_ref(copy.deepcopy(ref), bin_input[i].copy()))


# plot histograms of difference between instrument and manual bin        
fig4, axs4 = plt.subplots(1,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig4.subplots_adjust(hspace = 1, wspace=.001)
axs4 = axs4.ravel()

# start and stop for plotting
rstart, rstop = 1, -1
cstart,cstop = 1, -1

# bin (histogram) size
binn = 30

for i in range(0,len(CCDd_list)):
            
    ccdimg = CCDs_sub_img[i].copy()
    diff_img = ccdimg[rstart:rstop,cstart:cstop].copy()-binned[i][rstart:rstop,cstart:cstop].copy()
    
    meand = diff_img.mean()
    stdd = diff_img.std()
    #range=(meand-3*stdd,meand+3*stdd)

    axs4[i].hist(diff_img.ravel(), bins=binn, alpha=0.6,
                color = "skyblue", label='manual', density=True)
    
    axs4[i].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns'])+' man')
    
    fig4.suptitle('histogram of instrument-manual')
    

# plot histograms of manual and instrument bins
fig3, axs3 = plt.subplots(3,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig3.subplots_adjust(hspace = 1, wspace=.001)
axs3 = axs3.ravel()

# start and stop for plotting
rstart, rstop = 1, -1
cstart,cstop = 1, -1

# bin (histogram) size
binn = 80

for i in range(0,len(CCDd_list)):
        
    ccdimg = CCDs_sub_img[i].copy()
    
    if list_name == 'column.txt':
        fig3.suptitle('manual binning of shutter images (dark image subtracted) (on-chip column bin)')
    
    if list_name == 'row.txt':
        fig3.suptitle('manual binning of shutter images (dark image subtracted)')
        
    if list_name == 'fpga.txt':
        fig3.suptitle('manual binning of shutter images (dark image subtracted) (FPGA)')
        

    meanb = binned[i].mean()
    stdb = binned[i].std()
    
    meanc = ccdimg.mean()
    stdc = ccdimg.std()
    
    

    axs3[i].hist(binned[i][rstart:rstop,cstart:cstop].ravel(),range=(meanb-3*stdb,meanb+3*stdb), bins=binn, alpha=0.6,
                color = "skyblue", label='manual', density=True)
    axs3[i+len(CCDd_list)].hist(ccdimg[rstart:rstop,cstart:cstop].ravel(), range=(meanc-3*stdc,meanc+3*stdc), bins=binn, alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDd_list)].hist(ccdimg[rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(meanb-8*stdb,meanb+8*stdb), alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDd_list)].hist(binned[i][rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(meanb-8*stdb,meanb+8*stdb), alpha=0.6, 
                color='skyblue', label='manual', density=True) 
    axs3[i].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns'])+' man')
    axs3[i+len(CCDd_list)].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns']) + ' inst')
    axs3[i+2*len(CCDd_list)].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NCBIN CCDColumns']*CCDd_list[i]['NCBIN FPGAColumns']) + ' comp')




