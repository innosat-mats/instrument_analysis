#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 12:23:19 2020

@author: bjorn
"""



import sys
sys.path.append('/home/olemar/Projects/MATS/MATS-L0-processing')
sys.path.append('/home/olemar/Projects/MATS/instrument_analysis/read_images')

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
    
        # simple code for binning
        
        nrow, ncol, nrskip, ncskip, nrbin, ncbin = (ccd['NROW'], ccd['NCOL']+1,
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], ccd['NColBinCCD'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], ref['NColBinCCD'])
    
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
                ccd['NRSKIP'], ccd['NCSKIP'], ccd['NRBIN'], 2**ccd['NColBinFPGA'])
        
        nrowr, ncolr, nrskipr, ncskipr, nrbinr, ncbinr = (ref['NROW'], ref['NCOL']+1,
                ref['NRSKIP'], ref['NCSKIP'], ref['NRBIN'], 2**ref['NColBinFPGA'])
    
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
    
    return image1.copy()-image2.copy()



cal_day = '/080720_binning_nonlinearity'

# this is the reference exposure time (s)
ref_time = 4

# list containing file names of exposure time tests
list_exposure = 'exposure.txt'



####################################################
############# Load Reference Images ################
####################################################

dirname = ('/home/olemar/Projects/MATS/MATS-data/PayloadImages_200807_2')

os.chdir(dirname)

meanb_row, meanb_column, meanb_fpga = [],[],[]
meanc_row, meanc_column, meanc_fpga = [],[],[]
bins_row, bins_column, bins_fpga = [], [], []
means_row, means_column, means_fpga, means_exp = [],[],[],[]


# lists containing file names of binning tests
list_names = ['row.txt', 'column.txt', 'fpga.txt']



for list_name in list_names:
    
    CCDitems = []
    binned = []
    
    itemlist = np.genfromtxt(list_name, dtype='str')
           
    for line in itemlist:
        
        CCD = read_CCDitem_from_imgview(dirname, line)
            
        CCDitems.append(copy.deepcopy(CCD))
    
    CCDd_list = copy.deepcopy(CCDitems[0::2])
    CCDs_list = copy.deepcopy(CCDitems[1::2])
    
    
    # SUBTRACT DARK IMAGES
    CCDs_sub_img = []
    
    for i in range(0,len(CCDd_list)):
        CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))

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
            binned.append(bin_ref_FPGA(copy.deepcopy(ref), copy.deepcopy(bin_input[i])))
        
        else:
            binned.append(bin_ref(copy.deepcopy(ref), copy.deepcopy(bin_input[i])))
    
    
    # calculate means and arrays for plotting
    for i in range(0,len(CCDd_list)):
        
        ccdimg = CCDs_sub_img[i].copy()
        
        if list_name == 'fpga.txt':
            meanb_fpga.append(binned[i].mean())
            meanc_fpga.append(ccdimg.mean())
            bins_fpga.append(2**bin_input[i]['NColBinFPGA'])
            
            # with dark
            means_fpga.append(CCDs_list[i]['IMAGE'].mean())
            

        if list_name == 'row.txt':
            meanb_row.append(binned[i].mean())
            meanc_row.append(ccdimg.mean())
            bins_row.append(bin_input[i]['NRBIN'])
            
            # with dark
            means_row.append(CCDs_list[i]['IMAGE'].mean())

            
        if list_name == 'column.txt':
            meanb_column.append(binned[i].mean())
            meanc_column.append(ccdimg.mean())
            bins_column.append(bin_input[i]['NColBinCCD'])
            
            # with dark
            means_column.append(CCDs_list[i]['IMAGE'].mean())


# convert to numpy arrays
meanb_fpga = np.array(meanb_fpga)
meanc_fpga = np.array(meanc_fpga)

meanb_column = np.array(meanb_column)
meanc_column = np.array(meanc_column)

meanb_row = np.array(meanb_row)
meanc_row = np.array(meanc_row)
        

# EXPOSURE TIME TESTS

CCDitems = []
itemlist = np.genfromtxt(list_exposure, dtype='str')

for line in itemlist:
    CCD = read_CCDitem_from_imgview(dirname, line)
    CCDitems.append(copy.deepcopy(CCD))

CCDd_list = copy.deepcopy(CCDitems[0::2])
CCDs_list = copy.deepcopy(CCDitems[1::2])

# SUBTRACT DARK IMAGES

CCDs_sub_img, exp_times = [], []
mean_exp = []

for i in range(0,len(CCDd_list)):

    CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))
    exp_times.append(CCDs_list[i]['TEXPMS'])
    mean_exp.append(CCDs_sub_img[i].mean())
    
    # with dark
    means_exp.append(CCDs_list[i]['IMAGE'].mean())

# convert to numpy arrays    
exp_times = np.array(exp_times)
mean_exp = np.array(mean_exp)

# convert from TEXPMS to actual exposure time (s)
exp_times = ((exp_times - 2000)/1000)
exp_times = exp_times.astype(int)

# estimate linear response to increased exposure time
mult_exp = CCDs_sub_img[1].mean() * exp_times.copy()/4



# for removing any saturated measurements

bins_fpga = bins_fpga[0:-1]
meanc_fpga = meanc_fpga[0:-1]
meanb_fpga = meanb_fpga[0:-1]
means_fpga = means_fpga[0:-1]

bins_column = bins_column[0:-2]
meanc_column = meanc_column[0:-2]
meanb_column = meanb_column[0:-2]
means_column = means_column[0:-2]

bins_row = bins_row[0:-1]
meanc_row = meanc_row[0:-1]
meanb_row = meanb_row[0:-1]
means_row = means_row[0:-1]



####################################################
############# PLOTTING AND FITTING  ################
####################################################

# marker size
msize = 12

plt.figure(1)    


coef = np.polyfit(bins_fpga,meanc_fpga,1)
print(coef)
poly1d_fn = np.poly1d(coef) 
plt.plot(bins_fpga, poly1d_fn(bins_fpga), linewidth=1, color='blue', linestyle='--')
#plt.scatter(bins_fpga, meanb_fpga, color='blue', s=msize, label='manual fpga bin')
plt.scatter(bins_fpga, meanc_fpga, color='blue',s=msize, marker='x', label='fpga')

coef = np.polyfit(bins_row,meanc_row,1)
poly1d_fn = np.poly1d(coef) 
print(coef)
plt.plot(bins_row, poly1d_fn(bins_row), linewidth=1, color='red', linestyle='--')
#plt.scatter(bins_row, meanb_row, color='red',s=msize, label='manual row bin')
plt.scatter(bins_row, meanc_row, color='red',s=msize, marker='x', label='on-chip row bin')

coef = np.polyfit(bins_column,meanc_column,1)
poly1d_fn = np.poly1d(coef)
print(coef) 
plt.plot(bins_column, poly1d_fn(bins_column), linewidth=1, color='green', linestyle='--')
plt.scatter(bins_column, meanb_column, color='green',s=msize, label='manual on-chip column bin')
plt.scatter(bins_column, meanc_column, color='green', s=msize,marker='x', label='on-chip column bin')

# coef = np.polyfit(exp_times/ref_time,mean_exp,1)
# poly1d_fn = np.poly1d(coef) 
# print(coef)
# plt.plot(exp_times/ref_time, poly1d_fn(exp_times/ref_time), linewidth=1, color='pink', linestyle='--')
# plt.scatter(exp_times/ref_time, mean_exp, color='pink', s=msize, marker='x', label='int')
# plt.scatter(exp_times/ref_time, mult_exp, color='pink', s=msize,label='linear int')

#plt.xlabel('factor')
plt.title('mean pixel value (IR1)')
plt.legend()
plt.grid()
            


plt.figure(2)

plt.scatter(bins_fpga, meanc_fpga/meanb_fpga, color='blue', label='fpga')
plt.scatter(bins_row, meanc_row/meanb_row, color='red', label='row')
plt.scatter(bins_column, meanc_column/meanb_column, color='green', label='column')
plt.scatter(exp_times/ref_time, mean_exp/mult_exp, color='pink',  marker='x', label='int')

plt.ylim([0.965,1.005])
plt.xlim([0,70])

plt.title('mean pixel value measured/expected (IR1)')
#plt.xlabel('bin size (or exposure time factor)')

plt.legend()
    
plt.grid(2)


# with dark currnt
plt.figure(3)    


coef = np.polyfit(bins_fpga,means_fpga,1)
poly1d_fn = np.poly1d(coef) 
plt.plot(bins_fpga, poly1d_fn(bins_fpga), linewidth=1, color='blue', linestyle='--')
#plt.scatter(bins_fpga, meanb_fpga, color='blue', s=msize, label='manual fpga bin')
plt.scatter(bins_fpga, means_fpga, color='blue',s=msize, marker='x', label='fpga')

coef = np.polyfit(bins_row,means_row,1)
poly1d_fn = np.poly1d(coef) 
plt.plot(bins_row, poly1d_fn(bins_row), linewidth=1, color='red', linestyle='--')
#plt.scatter(bins_row, meanb_row, color='red',s=msize, label='manual row bin')
plt.scatter(bins_row, means_row, color='red',s=msize, marker='x', label='on-chip row bin')

coef = np.polyfit(bins_column,means_column,1)
poly1d_fn = np.poly1d(coef) 
plt.plot(bins_column, poly1d_fn(bins_column), linewidth=1, color='green', linestyle='--')
#plt.scatter(bins_column, meanb_column, color='green',s=msize, label='manual on-chip column bin')
plt.scatter(bins_column, means_column, color='green', s=msize,marker='x', label='on-chip column bin')

coef = np.polyfit(exp_times/4,means_exp,1)
poly1d_fn = np.poly1d(coef) 
plt.plot(exp_times/4, poly1d_fn(exp_times/4), linewidth=1, color='pink', linestyle='--')
plt.scatter(exp_times/4, means_exp, color='pink', s=msize, marker='x', label='int')
#plt.scatter(exp_times/4, mult_exp, color='pink', s=msize,label='linear int')

#plt.xlabel('factor')
plt.title('mean pixel value (inc dark) (IR1)')
plt.legend()
plt.grid()
            
plt.show(3)
    

    
    


