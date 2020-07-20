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
from LindasCalibrationFunctions import plotCCDitem
from LindasCalibrationFunctions import plot_simple
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics
from L1_calibration_functions import get_true_image_remove, desmear_true_image_remove, desmear_true_image
import copy

def bin_ref(ref, ccd):
    
        ## simple code for binning, does not take into consideration images where
        ## bad columns or rows have been omitted.
        nrow, ncol, nrskip, ncskip, nrbin, ncbin = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NColBinCCD'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NColBinCCD'])
    
        # reference image that will be binned according to ccd
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
    
        ## simple code for binning, does not take into consideration images where
        ## bad columns or rows have been omitted.
        nrow, ncol, nrskip, ncskip, nrbin, ncbin = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], 2**ccd['NColBinFPGA'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], 2**ref['NColBinFPGA'])
    
        # reference image that will be binned according to ccd
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



cal_day = '/080720_binning_nonlinearity'
#cal_day = '/20200713BinningTest'

####################################################
############# Load Reference Images ################
####################################################

dirname = ('/Users/bjorn/Documents/PhD/MATS/calibration'+ cal_day 
                + '/PayloadImages')

os.chdir(dirname)

CCDitems = []
binned = []

#list_names = ['row.txt', 'column.txt', 'fpga.txt']

itemlist = np.genfromtxt('row.txt', dtype='str')
#itemlist = np.genfromtxt('column.txt', dtype='str')
#itemlist = np.genfromtxt('column_11.txt', dtype='str')
#itemlist = np.genfromtxt('exposure.txt', dtype='str')
#itemlist = np.genfromtxt('fpga.txt', dtype='str')


       
for line in itemlist:
    CCD = read_CCDitem_from_imgview(dirname, line)
        
    CCDitems.append(CCD)

label = 'shutter'
# k = 0 on chip cbin, k=12 row bin, k=24 fpga
k=12

CCDd_list = np.copy(CCDitems[0::2])
CCDs_list = np.copy(CCDitems[1::2])

# plot the obtained images

fig, axs = plt.subplots(2,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = 1, wspace=.001)
axs = axs.ravel()
fig.suptitle('images')
  
# LOOP TO PLOT IMAGES
for i in range(0,len(CCDd_list)):
    im = axs[i].imshow(CCDd_list[i]['IMAGE'],cmap='jet')
    mean=CCDd_list[i]['IMAGE'].mean()
    std=CCDd_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i].set_title('dark: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA']))
    
    im = axs[i+len(CCDd_list)].imshow(CCDs_list[i]['IMAGE'],cmap='jet')
    mean=CCDs_list[i]['IMAGE'].mean()
    std=CCDs_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i+len(CCDd_list)])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i+len(CCDd_list)].set_title('shutter: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA']))

# SUBTRACT DARK IMAGES
CCDs_sub_img = []

for i in range(0,len(CCDd_list)):

    #CCDns_sub_img.append(img_diff(CCDns_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))
    CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))
    
fig1, axs1 = plt.subplots(1,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig1.subplots_adjust(hspace = 1, wspace=.001)
axs1 = axs1.ravel()
fig1.suptitle('images with subtracted dark images')

# LOOP TO PLOT IMAGES
for i in range(0,len(CCDd_list)):
    im = axs1[i].imshow(CCDs_sub_img[i],cmap='jet')
    mean=CCDs_sub_img[i].mean()
    std=CCDs_sub_img[i].std()
    fig1.colorbar(im, ax=axs1[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs1[i].set_title('shutter - dark: ' + str(CCDd_list[i]['NRBIN']) + 'x' 
                      + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA']))
    


# BINNING

if label == 'shutter':
    bin_input = copy.deepcopy(CCDs_list)
    ref = copy.deepcopy(CCDs_list[0])
    ref['IMAGE'] = CCDs_sub_img[0].copy()

    for i in range(0,len(CCDd_list)): 
    
        bin_input[i]['IMAGE'] = CCDs_sub_img[i].copy()
    


for i in range(0,len(CCDd_list)): 
    
    if k == 24:
        binned.append(bin_ref_FPGA(copy.deepcopy(ref), bin_input[i].copy()))
    
    else:
        binned.append(bin_ref(copy.deepcopy(ref), bin_input[i].copy()))

        

fig4, axs4 = plt.subplots(1,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig4.subplots_adjust(hspace = 1, wspace=.001)
axs4 = axs4.ravel()

rstart, rstop = 1, -1
cstart,cstop = 1, -1

for i in range(0,len(CCDd_list)):
            
    
    if label == 'shutter':
        ccdimg = CCDs_sub_img[i].copy()
        
    binn = 30
    
    diff_img = ccdimg[rstart:rstop,cstart:cstop].copy()-binned[i][rstart:rstop,cstart:cstop].copy()
    
    meand = diff_img.mean()
    stdd = diff_img.std()
    #range=(meand-3*stdd,meand+3*stdd)

    axs4[i].hist(diff_img.ravel(), bins=binn, alpha=0.6,
                color = "skyblue", label='manual', density=True)
    
    axs4[i].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA'])+' man')
    
    fig4.suptitle('histogram of instrument-manual')
    

fig3, axs3 = plt.subplots(3,len(CCDd_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig3.subplots_adjust(hspace = 1, wspace=.001)
axs3 = axs3.ravel()

rstart, rstop = 1, -1
cstart,cstop = 1, -1

for i in range(0,len(CCDd_list)):
            
    
    if label == 'shutter':
        ccdimg = CCDs_sub_img[i].copy()
        
        if k == 0:
            fig3.suptitle('manual binning of shutter images (dark image subtracted) (on-chip column bin)')
        
        if k == 12:
            fig3.suptitle('manual binning of shutter images (dark image subtracted)')
            
        if k == 24:
            fig3.suptitle('manual binning of shutter images (dark image subtracted) (FPGA)')
        
    
        
    #print(ccdimg)
    meanb = binned[i].mean()
    stdb = binned[i].std()
    
    meanc = ccdimg.mean()
    stdc = ccdimg.std()
    
    binn = 80

    axs3[i].hist(binned[i][rstart:rstop,cstart:cstop].ravel(),range=(meanb-3*stdb,meanb+3*stdb), bins=binn, alpha=0.6,
                color = "skyblue", label='manual', density=True)
    axs3[i+len(CCDd_list)].hist(ccdimg[rstart:rstop,cstart:cstop].ravel(), range=(meanc-3*stdc,meanc+3*stdc), bins=binn, alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDd_list)].hist(ccdimg[rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(meanb-8*stdb,meanb+8*stdb), alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDd_list)].hist(binned[i][rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(meanb-8*stdb,meanb+8*stdb), alpha=0.6, 
                color='skyblue', label='manual', density=True) 
    axs3[i].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA'])+' man')
    axs3[i+len(CCDd_list)].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA']) + ' inst')
    axs3[i+2*len(CCDd_list)].set_title(str(CCDd_list[i]['NRBIN']) + 'x'
                        + str(CCDd_list[i]['NColBinCCD']*2**CCDd_list[i]['NColBinFPGA']) + ' comp')




