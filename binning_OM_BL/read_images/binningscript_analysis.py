#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:23:19 2020

@author: bjorn
"""



import sys
sys.path.append('/home/olemar/Projects/MATS/MATS-L0-processing')
sys.path.append('/home/olemar/Projects/MATS/MATS-L1-processing')

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from LindasCalibrationFunctions import plotCCDitem
from LindasCalibrationFunctions import plot_simple
from read_in_functions import read_CCDitem_from_imgview, readimageviewpics
from L1_calibration_functions import get_true_image_remove, desmear_true_image_remove, desmear_true_image
import copy


# TO DO: COMBINE THE FOLLOWING TWO FUNCITONS INTO ONE THAT BINS ACCORDING TO
# BOTH FPGA AND ON-CHIP & ROW SETTINGS; 

def bin_ref(ref, ccd):
    
        # simple code for binning
        
        nrow, ncol, nrskip, ncskip, nrbin, ncbin, exptime = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NCBIN CCDColumns'],ccd['TEXPMS'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr, exptimer = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NCBIN CCDColumns'],ref['TEXPMS'])
    
        exptimefactor = int((exptime-2000)/(exptimer-2000))
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
                
            
            binned = binned*exptimefactor
            return binned

        else:
            
            sys.exit('Error: images not from the same CCD region.')




def bin_ref_FPGA(ref, ccd):
    
        # simple code for binning 
        nrow, ncol, nrskip, ncskip, nrbin, ncbin, exptime = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NCBIN CCDColumns'],ccd['TEXPMS'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr, exptimer = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NCBIN CCDColumns'],ref['TEXPMS'])
    
        exptimefactor = int((exptime-2000)/(exptimer-2000))
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
           
            binned = binned*exptimefactor
            return binned

        else:
            
            sys.exit('Error: images not from the same CCD region.')
            
            
def img_diff(image1, image2):
    
    return image1-image2



####################################################
#############        LOAD DATA      ################
####################################################


cal_day = '20200812'
# type of binning 
channel = 5
#binning_type = 'column'
#binning_type = 'row'
binning_type = 'exptime'

dirname = ('/home/olemar/Projects/MATS/MATS-data/binning_test_channel_'  + str(channel) + '_' + cal_day + '/' + binning_type)

CCDitems=[]
IDstrings=[]
binned=[]


os.chdir(dirname)

for file in glob.glob("*.pnm"):
    IDstrings.append(file.strip('.pnm'))
    
# sort by time
IDstrings.sort(key=int)

for IDstring in IDstrings:    
    CCD = read_CCDitem_from_imgview(dirname, IDstring)
    CCDitems.append(CCD)    


# long exposure, short exposure, new reference images
CCDl_list = np.copy(CCDitems[0::4])
CCDs_list = np.copy(CCDitems[1::4])
CCDr_list = np.copy(CCDitems[2::4])
CCDrs_list = np.copy(CCDitems[3::4]) #reference short (never binned)

####################################################
#############       PLOTTING 1   ##################
####################################################


# plot the images
fig, axs = plt.subplots(2,len(CCDs_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = 1, wspace=.001)
axs = axs.ravel()
fig.suptitle('images')
  
for i in range(0,len(CCDs_list)):
    im = axs[i].imshow(CCDs_list[i]['IMAGE'],cmap='jet')
    mean=CCDs_list[i]['IMAGE'].mean()
    std=CCDs_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i].set_title('short: ' + str(CCDs_list[i]['NRBIN']) + 'x' 
                      + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns']))
    
    im = axs[i+len(CCDs_list)].imshow(CCDl_list[i]['IMAGE'],cmap='jet')
    mean=CCDl_list[i]['IMAGE'].mean()
    std=CCDl_list[i]['IMAGE'].std()
    fig.colorbar(im, ax=axs[i+len(CCDs_list)])
    im.set_clim(mean-2*std,mean+2*std)
    axs[i+len(CCDs_list)].set_title('long: ' + str(CCDs_list[i]['NRBIN']) + 'x' 
                      + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns']))

# SUBTRACT DARK IMAGES
CCDl_sub_img,CCDr_sub_img = [],[]

for i in range(0,len(CCDs_list)):
    
    print(i)
    
    # subtract dark current from both long and references
    CCDl_sub_img.append(img_diff(CCDl_list[i]['IMAGE'].copy(), CCDs_list[i]['IMAGE'].copy()))
    CCDr_sub_img.append(img_diff(CCDr_list[i]['IMAGE'].copy(), CCDrs_list[i]['IMAGE'].copy())) # update

# plot the images with removed dark current    
fig1, axs1 = plt.subplots(1,len(CCDs_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig1.subplots_adjust(hspace = 1, wspace=.001)
axs1 = axs1.ravel()
fig1.suptitle('images with subtracted dark images')

for i in range(0,len(CCDs_list)):
    im = axs1[i].imshow(CCDl_sub_img[i],cmap='jet')
    mean=CCDl_sub_img[i].mean()
    std=CCDl_sub_img[i].std()
    fig1.colorbar(im, ax=axs1[i])
    im.set_clim(mean-2*std,mean+2*std)
    axs1[i].set_title('long - short: ' + str(CCDs_list[i]['NRBIN']) + 'x' 
                      + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns']))
    

####################################################
#############        BINNING       #################
####################################################


# copy settings from long image (for binning settings)
bin_input = copy.deepcopy(CCDl_list)


# replace images with the images with subtrated dark (I think this is old and not used)
for i in range(0,len(CCDs_list)): 
    bin_input[i]['IMAGE'] = CCDl_sub_img[i].copy()
    

# create manually binned images
for i in range(0,len(CCDs_list)): 
    
    # update reference image that should be binned manually
    ref = copy.deepcopy(CCDr_list[i])
    ref['IMAGE'] = CCDr_sub_img[i].copy()
    
    
    # bin reference image according to bin_input settings
    if binning_type == 'fpga':
        binned.append(bin_ref_FPGA(copy.deepcopy(ref), bin_input[i].copy()))
    
    else:
        binned.append(bin_ref(copy.deepcopy(ref), bin_input[i].copy()))



####################################################
#############       PLOTTING 2    ##################
####################################################

     
# fig4, axs4 = plt.subplots(1,len(CCDs_list), figsize=(15, 6), facecolor='w', edgecolor='k')
# fig4.subplots_adjust(hspace = 1, wspace=.001)
# axs4 = axs4.ravel()

# # start and stop for plotting
# rstart, rstop = 1, -1
# cstart,cstop = 1, -1

# # histogram number of bins
# binn = 200

# for i in range(0,len(CCDs_list)):
            
#     # # plot 2D difference images
#     # inst_bin = CCDl_sub_img[i].copy()
#     # diff_img = inst_bin[rstart:rstop,cstart:cstop].copy()-binned[i][rstart:rstop,cstart:cstop].copy()
#     # mean=diff_img.mean()
#     # std=diff_img.std()
#     # im = axs4[i].imshow(diff_img,cmap='jet')
#     # fig4.colorbar(im, ax=axs4[i])
#     # im.set_clim(mean-2*std,mean+2*std)
#     # #im.set_clim(0.97,1.03)

#     # fig4.colorbar(im, ax=axs[i+len(CCDs_list)])
#     # fig4.suptitle('instrument-manual')


#     # plot histograms of difference between instrument and manual bin   
#     inst_bin = CCDl_sub_img[i].copy()
#     diff_img = inst_bin[rstart:rstop,cstart:cstop].copy()-binned[i][rstart:rstop,cstart:cstop].copy()
#     mean = diff_img.mean()
#     std = diff_img.std()
#     axs4[i].hist(diff_img.ravel(), bins=binn, alpha=0.6,
#                 color = "skyblue", label='manual', density=True, range=(mean-3*std,mean+3*std))
#     axs4[i].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
#                         + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns'])+' man')
#     fig4.suptitle('histogram of instrument-manual')
    

# plot histograms of manual and instrument bins
fig3, axs3 = plt.subplots(3,len(CCDs_list), figsize=(15, 6), facecolor='w', edgecolor='k')
fig3.subplots_adjust(hspace = 1, wspace=.001)
axs3 = axs3.ravel()

# start and stop for plotting
rstart, rstop = 1, -1
cstart,cstop = 1, -1

# bin (histogram) size
binn = 80

for i in range(0,len(CCDs_list)):
        
    inst_bin = CCDl_sub_img[i].copy()
    
    if binning_type == 'column':
        fig3.suptitle('manual binning (dark img subtracted) (on-chip column bin)')
    
    if binning_type == 'row':
        fig3.suptitle('manual binning (dark img subtracted)')
        
    if binning_type == 'fpga':
        fig3.suptitle('manual binning (dark img subtracted) (FPGA)')
        

    mean_man = binned[i].mean()
    std_man = binned[i].std()
    
    mean_inst = inst_bin.mean()
    std_inst = inst_bin.std()
    
    axs3[i].hist(binned[i][rstart:rstop,cstart:cstop].ravel(),range=(mean_man-3*std_man,mean_man+3*std_man), bins=binn, alpha=0.6,
                color = "skyblue", label='manual', density=True)
    axs3[i+len(CCDs_list)].hist(inst_bin[rstart:rstop,cstart:cstop].ravel(), range=(mean_inst-3*std_man,mean_inst+3*std_inst), bins=binn, alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDs_list)].hist(inst_bin[rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(mean_man-3*std_man,mean_man+3*std_man), alpha=0.4, 
                color='red', label='instrument', density=True) 
    axs3[i+2*len(CCDs_list)].hist(binned[i][rstart:rstop,cstart:cstop].ravel(), bins=binn, range=(mean_man-3*std_man,mean_man+3*std_man), alpha=0.6, 
                color='skyblue', label='manual', density=True) 
    axs3[i].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
                        + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns'])+' man')
    axs3[i+len(CCDs_list)].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
                        + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns']) + ' inst')
    axs3[i+2*len(CCDs_list)].set_title(str(CCDs_list[i]['NRBIN']) + 'x'
                        + str(CCDs_list[i]['NCBIN CCDColumns']*CCDs_list[i]['NCBIN FPGAColumns']) + ' comp')




