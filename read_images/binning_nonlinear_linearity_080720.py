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
    
    return image1.copy()-image2.copy()



cal_day = '/080720_binning_nonlinearity'
#cal_day = '/20200713BinningTest'

ref_time = 4
list_exposure = 'exposure.txt'
#list_exposure = 'exposure_11.txt'


####################################################
############# Load Reference Images ################
####################################################

dirname = ('/Users/bjorn/Documents/PhD/MATS/calibration'+ cal_day 
                + '/PayloadImages')

os.chdir(dirname)

meanb_row, meanb_column, meanb_fpga = [],[],[]
meanc_row, meanc_column, meanc_fpga = [],[],[]
bins_row, bins_column, bins_fpga = [], [], []
means_row, means_column, means_fpga, means_exp = [],[],[],[]

list_names = ['row.txt', 'column.txt', 'fpga.txt']
#list_names = ['row_11.txt', 'column_11.txt', 'fpga.txt']
#list_names = ['column.txt']

for list_name in list_names:
    
    CCDitems = []
    binned = []
    
    itemlist = np.genfromtxt(list_name, dtype='str')

    if list_name == 'fpga.txt':
        k = 24
    else:
        k = 0
           
    for line in itemlist:
        
        CCD = read_CCDitem_from_imgview(dirname, line)
            
        CCDitems.append(copy.deepcopy(CCD))
    
    label = 'shutter'
    
    CCDd_list = copy.deepcopy(CCDitems[0::2])
    CCDs_list = copy.deepcopy(CCDitems[1::2])
    
    
    # SUBTRACT DARK IMAGES
    CCDs_sub_img = []
    
    for i in range(0,len(CCDd_list)):
    
        CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))

    
    # BINNING
    
    if label == 'shutter':
        bin_input = copy.deepcopy(CCDs_list)
        ref = copy.deepcopy(CCDs_list[0])
        ref['IMAGE'] = CCDs_sub_img[0].copy()
    
        for i in range(0,len(CCDd_list)): 
        
            bin_input[i]['IMAGE'] = CCDs_sub_img[i].copy()
        
    
    
    for i in range(0,len(CCDd_list)): 
        
        if k == 24:
            binned.append(bin_ref_FPGA(copy.deepcopy(ref), copy.deepcopy(bin_input[i])))
        
        else:
            binned.append(bin_ref(copy.deepcopy(ref), copy.deepcopy(bin_input[i])))
    
    
    for i in range(0,len(CCDd_list)):
        
        
        if label == 'shutter':
            ccdimg = CCDs_sub_img[i].copy()

        #print(ccdimg)
        
        #stdb = binned[i].std()
        #stdc = ccdimg.std()
        

        
        if list_name == 'fpga.txt':
            meanb_fpga.append(binned[i].mean())
            meanc_fpga.append(ccdimg.mean())
            bins_fpga.append(2**bin_input[i]['NColBinFPGA'])
            
            means_fpga.append(CCDs_list[i]['IMAGE'].mean())
            

        if list_name == 'row.txt':
            meanb_row.append(binned[i].mean())
            meanc_row.append(ccdimg.mean())
            bins_row.append(bin_input[i]['NRBIN'])
            
            means_row.append(CCDs_list[i]['IMAGE'].mean())

            
        if list_name == 'column.txt':
            meanb_column.append(binned[i].mean())
            meanc_column.append(ccdimg.mean())
            bins_column.append(bin_input[i]['NColBinCCD'])
            
            means_column.append(CCDs_list[i]['IMAGE'].mean())


meanb_fpga = np.array(meanb_fpga)
meanc_fpga = np.array(meanc_fpga)

meanb_column = np.array(meanb_column)
meanc_column = np.array(meanc_column)

meanb_row = np.array(meanb_row)
meanc_row = np.array(meanc_row)
        

### LOAD EXP IMAGES

CCDitems = []

itemlist = np.genfromtxt(list_exposure, dtype='str')

       
for line in itemlist:
    
    CCD = read_CCDitem_from_imgview(dirname, line)
        
    CCDitems.append(copy.deepcopy(CCD))

label = 'shutter'

CCDd_list = copy.deepcopy(CCDitems[0::2])
CCDs_list = copy.deepcopy(CCDitems[1::2])


# SUBTRACT DARK IMAGES & ADD EXP TIMES
CCDs_sub_img, exp_times = [], []
mean_exp = []

for i in range(0,len(CCDd_list)):

    CCDs_sub_img.append(img_diff(CCDs_list[i]['IMAGE'].copy(), CCDd_list[i]['IMAGE'].copy()))
    exp_times.append(CCDs_list[i]['TEXPMS'])
    mean_exp.append(CCDs_sub_img[i].mean())
    
    means_exp.append(CCDs_list[i]['IMAGE'].mean())
    
exp_times = np.array(exp_times)
mean_exp = np.array(mean_exp)
exp_times = ((exp_times - 2000)/1000)
exp_times = exp_times.astype(int)

mult_exp = CCDs_sub_img[1].mean() * exp_times.copy()/4



# coef = np.polyfit(bins_fpga,meanb_fpga,1)
# poly1d_fn = np.poly1d(coef) 

# coef = np.polyfit(bins_fpga,meanc_fpga,1)
# poly1d_fn = np.poly1d(coef) 

# coef = np.polyfit(bins_row,meanb_row,1)
# poly1d_fn = np.poly1d(coef) 

# coef = np.polyfit(bins_row,meanc_row,1)
# poly1d_fn = np.poly1d(coef) 



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
    

    
    


