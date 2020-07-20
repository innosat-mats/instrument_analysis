#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 17:09:57 2019

@author: franzk, Linda Megner
(Linda has removed some fuctions and added more . Original franzk script is L1_functions.py)


Functions used for MATS L1 processing, based on corresponding MATLAB scripts provided by
Georgi Olentsenko and Mykola Ivchenko
The MATLAB script can be found here: https://github.com/OleMartinChristensen/MATS-image-analysis



"""

import numpy as np
import scipy.io

# some_file.py
#import sys
# insert at 1, 0 is the script path (or '' in REPL)


#sys.path.insert(1, '/Users/lindamegner/MATS/retrieval/Level0/MATS-L0-processing-master')

#ismember function taken from https://stackoverflow.com/questions/15864082/python-equivalent-of-matlabs-ismember-function
#needs to be tested, if actually replicates matlab ismember function


class CCD:
    def __init__(self, channel):
        self.channel=channel
        if channel=='IR1':
            CCDID=16
        elif channel=='IR2':
            CCDID=17
        elif channel=='IR3':
            CCDID=18
        elif channel=='IR4':
            CCDID=19
        elif channel=='NADIR':
            CCDID=20
        elif channel=='UV1':
            CCDID=21
        elif channel=='UV2':
            CCDID=22            
        elif channel=='KTH test channel':
            CCDID=16 
        

        filename='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_KTH/MATS_CCD_DC_calibration_FINAL/FM0'+str(CCDID)+'_CCD_DC_calibration.mat'

#        darkdir='/Users/lindamegner/MATS/retrieval/Calibration/FM_calibration_at_KTH/FM_CCD_DC_calibration_data_reduced/'
#        filename= darkdir+ 'FM0' +str(CCDID) +'_CCD_DC_calibration_DATA_reduced.mat'
        mat = scipy.io.loadmat(filename)

        self.dc_zero_avr_HSM=mat['dc_zero_avr_HSM']
        self.dc_zero_std_HSM=mat['dc_zero_std_HSM']
        self.dc_zero_avr_LSM=mat['dc_zero_avr_LSM']
        self.dc_zero_std_LSM=mat['dc_zero_std_LSM']
        
        self.image_HSM=mat['image_HSM']
        self.image_LSM=mat['image_LSM']
        
        
        self.ro_avr_HSM=mat['ro_avr_HSM']
        self.ro_std_HSM=mat['ro_std_HSM']
        self.alpha_avr_HSM=mat['alpha_avr_HSM']
        self.alpha_std_HSM=mat['alpha_std_HSM']
        
        self.ro_avr_LSM=mat['ro_avr_LSM']
        self.ro_std_LSM=mat['ro_std_LSM']
        self.alpha_avr_LSM=mat['alpha_avr_LSM']
        self.alpha_std_LSM=mat['alpha_std_LSM']
        
        self.log_a_avr_HSM=mat['log_a_avr_HSM']
        self.log_a_std_HSM=mat['log_a_std_HSM']
        self.log_b_avr_HSM=mat['log_b_avr_HSM']
        self.log_b_std_HSM=mat['log_b_std_HSM']
        
        self.log_a_avr_LSM=mat['log_a_avr_LSM']
        self.log_a_std_LSM=mat['log_a_std_LSM']
        self.log_b_avr_LSM=mat['log_b_avr_LSM']
        self.log_b_std_LSM=mat['log_b_std_LSM']
        
        #2D dark current subtraction stuff
        self.log_a_img_avr_LSM=mat['log_a_img_avr_LSM']
        self.log_a_img_err_LSM=mat['log_a_img_err_LSM']
        self.log_b_img_avr_LSM=mat['log_b_img_avr_LSM']
        self.log_b_img_err_LSM=mat['log_b_img_err_LSM']
        
        self.log_a_img_avr_HSM=mat['log_a_img_avr_HSM']
        self.log_a_img_err_HSM=mat['log_a_img_err_HSM']
        self.log_b_img_avr_HSM=mat['log_b_img_avr_HSM']
        self.log_b_img_err_HSM=mat['log_b_img_err_HSM']
        
        
        self.hot_pix=np.where(self.image_HSM>=0.8*np.max(self.image_HSM))

        if   (self.channel=='UV1' or self.channel=='UV2'):
            self.ampcorrection=3/2 #Amplification hack - check with Gabriel how to to properly
        else:
            self.ampcorrection=1
                  
    def darkcurrent(self, T, mode): #electrons/s

        if mode == 0: 
            darkcurrent=10**(self.log_a_avr_HSM*T+self.log_b_avr_HSM)
        elif mode == 1:
            darkcurrent=10**(self.log_a_avr_LSM*T+self.log_b_avr_LSM)
        else :
            print('Undefined mode')
        return darkcurrent
    
    def darkcurrent2D(self, T, mode): #electrons/s
        if mode == 0: 
            darkcurrent=10**(self.log_a_img_avr_HSM*T+self.log_b_img_avr_HSM)
        elif mode == 1:
            darkcurrent=10**(self.log_a_img_avr_LSM*T+self.log_b_img_avr_LSM)
        else :
            print('Undefined mode')
        return darkcurrent
    
    def ro_avr(self, mode):
        if mode == 0: 
            ro_avr=self.ro_avr_HSM
        elif mode == 1:
            ro_avr=self.ro_avr_LSM
        else :
            print('Undefined mode')
        return ro_avr
    
    def alpha_avr(self, mode): #electrons/LSB
        if mode == 0: 
            alpha_avr=self.alpha_avr_HSM
        elif mode == 1:
            alpha_avr=self.alpha_avr_LSM
        else :
            print('Undefined mode')
        return alpha_avr   

def subtract_dark(image, CCDitem):

    CCDunit=CCD(CCDitem['channel']) #TODO: Remember to clear out the non used bits of CCD /Linda
  
    totdarkcurrent=CCDunit.darkcurrent2D(CCDitem['temperature'],int(CCDitem['SigMode']))*int(CCDitem['TEXPMS'])/1000. # tot dark current in electrons
    totbinpix=int(CCDitem['NColBinCCD'])*2**int(CCDitem['NColBinFPGA']) #Note that the numbers are described in differnt ways see table 69 in Software iCD
    image_dark_sub=image-totbinpix*CCDunit.ampcorrection*totdarkcurrent/CCDunit.alpha_avr(int(CCDitem['SigMode']))
    return image_dark_sub


def winmode_correct(CCDitem):
    if (CCDitem['WinMode']) <=4:
        winfactor=2**CCDitem['WinMode']
    elif (CCDitem['WinMode']) ==7:       
        winfactor=1    
    else:
        raise Exception('Undefined Window')
    image_lsb=winfactor*CCDitem['IMAGE'] 
    return image_lsb

def get_true_image(image, header):
    #calculate true image by removing readout offset (pixel blank value) and
    #compensate for bad colums 
    
    ncolbinC=int(header['NColBinCCD'])
    if ncolbinC == 0:
        ncolbinC = 1
    
    #remove gain
    true_image = image * 2**(int(header['DigGain'])) 
    
    #bad column analysis
    n_read, n_coadd = binning_bc(int(header['NCOL'])+1, int(header['NCSKIP']), 2**int(header['NColBinFPGA']), ncolbinC, header['BC'])
    
    #go through the columns
    for j_c in range(0, int(header['NCOL'])):
        #remove blank values and readout offsets
        true_image[0:int(header['NROW']), j_c] = true_image[0:int(header['NROW']), j_c] - n_read[j_c]*(header['TBLNK']-128)-128
        
        #compensate for bad columns
        true_image[0:int(header['NROW']), j_c] = true_image[0:int(header['NROW']), j_c] * (2**int(header['NColBinFPGA'])*ncolbinC/n_coadd[j_c])
    
    return true_image


def get_true_image_remove(header):
    #calculate true image by removing readout offset (pixel blank value) and
    #compensate for bad colums 
    
    #note that header is now a CCDitem
    ncolbinC=int(header['NColBinCCD'])
    if ncolbinC == 0:
        ncolbinC = 1
    
    #remove gain
    true_image = header['IMAGE'] * 2**(int(header['DigGain'])) 
    
    #bad column analysis
    n_read, n_coadd = binning_bc(int(header['NCOL'])+1, int(header['NCSKIP']), 2**int(header['NColBinFPGA']), ncolbinC, header['BC'])

    #go through the columns
    for j_c in range(0, int(header['NCOL'])):
        #remove blank values and readout offsets
        true_image[0:int(header['NROW']), j_c] = true_image[0:int(header['NROW']), j_c] - n_read[j_c]*(header['TBLNK']-128)-128
        #true_image[0:int(header['NROW']), j_c] = true_image[0:int(header['NROW']), j_c] - n_read[j_c]*((header['TBLNK']+header['LBLNK'])/2-128)-128

        #compensate for bad columns
        true_image[0:int(header['NROW']), j_c] = true_image[0:int(header['NROW']), j_c] * (2**int(header['NColBinFPGA'])*ncolbinC/n_coadd[j_c])

    return true_image


def binning_bc(Ncol, Ncolskip, NcolbinFPGA, NcolbinCCD, BadColumns):
    
    #a routine to estimate the correction factors for column binning with bad columns
    
    #n_read - array, containing the number of individually read superpixels
    #           attributing to the given superpixel
    #n_coadd - array, containing the number of co-added individual pixels
    #Input - as per ICD. BadColumns - array containing the index of bad columns
    #           (the index of first column is 0)
    
    n_read=np.zeros(Ncol)
    n_coadd=np.zeros(Ncol)
    
    col_index=Ncolskip
    
    for j_col in range(0,Ncol):
        for j_FPGA in range(0,NcolbinFPGA):
            continuous=0
            for j_CCD in range(0,NcolbinCCD):
                if col_index in BadColumns:
                    if continuous == 1:
                        n_read[j_col]=n_read[j_col]+1
                    continuous=0
                else:
                    continuous=1
                    n_coadd[j_col]=n_coadd[j_col]+1
                
                col_index=col_index+1
            
            if continuous == 1:
                n_read[j_col]=n_read[j_col]+1
    
    return n_read, n_coadd


def desmear_true_image_remove(header):
    #note that header is now a CCDitem
    nrow = int(header['NROW'])
    ncol = int(header['NCOL']) + 1
    
    #calculate extra time per row
    T_row_extra, T_delay = calculate_time_per_row(header)
    T_exposure = int(header['TEXPMS'])/1000#check for results when shifting from python 2 to 3
    
    image=header['IMAGE']
    TotTime=0
    for irow in range(1,nrow):
        for krow in range(0,irow):
            image[irow,0:ncol]=image[irow,0:ncol] - image[krow,0:ncol]*(T_row_extra/T_exposure)
            TotTime=TotTime+T_row_extra
           
    #row 0 here is the first row to read out from the chip

    return image


def desmear_true_image(image, header):
    
    nrow = int(header['NROW'])
    ncol = int(header['NCOL']) + 1
    
    #calculate extra time per row
    T_row_extra, T_delay = calculate_time_per_row(header)
    T_exposure = int(header['TEXPMS'])/1000#check for results when shifting from python 2 to 3
    
    TotTime=0
    for irow in range(1,nrow):
        for krow in range(0,irow):
            image[irow,0:ncol]=image[irow,0:ncol] - image[krow,0:ncol]*(T_row_extra/T_exposure)
            TotTime=TotTime+T_row_extra
           
    #row 0 here is the first row to read out from the chip

    return image

def calculate_time_per_row(header):
    
    #this function provides some useful timing data for the CCD readout
    
    #Note that minor "transition" states may have been omitted resulting in 
    #somewhat shorter readout times (<0.1%).
    
    #Default timing setting is_
    #ccd_r_timing <= x"A4030206141D"&x"010303090313"
    
    #All pixel timing setting is the final count of a counter that starts at 0,
    #so the number of clock cycles exceeds the setting by 1
    
    #image parameters
    ncol=int(header['NCOL'])+1
    ncolbinC=int(header['NColBinCCD'])
    if ncolbinC == 0:
        ncolbinC = 1
    ncolbinF=2**int(header['NColBinFPGA'])
    
    nrow=int(header['NROW'])
    nrowbin=int(header['NRBIN'])
    if nrowbin == 0:
        nrowbin = 1
    nrowskip=int(header['NRSKIP'])
    
    n_flush=int(header['NFLUSH'])
    
    #timing settings
    #full_timing = 1 #TODO <-- meaning?
    full_timing = int(header['TimingFlag'])
    #full pixel readout timing n#TODO discuss this with OM,  LMc these are default values change these when the header contians this infromation
    
    time0 = 1 + 19 # x13%TODO
    time1 = 1 +  3 # x03%TODO
    time2 = 1 +  9 # x09%TODO
    time3 = 1 +  3 # x03%TODO
    time4 = 1 +  3 # x03%TODO
    time_ovl = 1 + 1 # x01%TODO
    
    # fast pixel readout timing
    timefast  = 1 + 2 # x02%TODO
    timefastr = 1 + 3 # x03%TODO
    
    #row shift timing
    row_step = 1 + 164 # xA4%TODO
    
    
    clock_period = 30.517 #master clock period, ns 32.768 MHz
    
    #there is one extra clock cycle, effectively adding to time 0
    Time_pixel_full = (1+ time0 + time1 + time2 + time3 + time4 + 3*time_ovl)*clock_period
    
    # this is the fast timing pixel period
    Time_pixel_fast = (1+ 4*timefast + 3*time_ovl + timefastr)*clock_period
    
    #here we calculate the number of fast and slow pixels
    #NOTE: the effect of bad pixels is disregarded here
    
    if full_timing == 1:
        n_pixels_full = 2148
        n_pixels_fast = 0
    else:
        if ncolbinC < 2: #no CCD binning
            n_pixels_full = ncol * ncolbinF
        else: #there are two "slow" pixels for one superpixel to be read out
            n_pixels_full = 2*ncol*ncolbinF
    
        n_pixels_fast = 2148 - n_pixels_full
    
    #time to read out one row
    T_row_read = n_pixels_full*Time_pixel_full + n_pixels_fast*Time_pixel_fast
    
    # shift time of a single row
    T_row_shift = (64 + row_step *10)*clock_period
    
    #time of the exposure start delay from the start_exp signal # n_flush=1023
    T_delay = T_row_shift * n_flush
    
    #total time of the readout
    T_readout = T_row_read*(nrow+nrowskip+1) + T_row_shift*(1+nrowbin*nrow)
    
    
    #"smearing time"
    #(this is the time that any pixel collects electrons in a wrong row, during the shifting.)
    #For smearing correction, this is the "extra exposure time" for each of the rows.
    
    T_row_extra = (T_row_read + T_row_shift*nrowbin) / 1e9    
    #T_row_extra = (T_row_read + T_row_shift*nrowbin) / 1e9    
    
    return T_row_extra, T_delay

def compare_image(image1, image2, header):
    
    #this is a function to compare two images of the same size
    #one comparison is a linear fit of columns, the other comparison is a linear fit
    #of rows, the third is a linear fit of the whole image
    
    sz1=image1.shape
    sz2=image2.shape
    
    if sz1[0] != sz2[0] or sz1[1] != sz2[1]:
        print('sizes of input images do not match')
    
    nrow=sz1[0]
    ncol=sz1[1]
    
    nrowskip = int(header['NRSKIP'])
    ncolskip = int(header['NCSKIP'])
    
    nrowbin = int(header['NRBIN'])
    ncolbinC = int(header['NCBIN'])
    ncolbinF = 2**int(header['NColBinFPGA'])
    
    if nrowskip + nrowbin*nrow > 511:
        nrow = np.floor((511-nrowskip)/nrowbin)
        
    if ncolskip + ncolbinC*ncolbinF*ncol > 2047:
        nrow = np.floor((2047-ncolskip)/(ncolbinC*ncolbinF))
    print(nrow,image1.shape)        
    image1 = image1[0:nrow-1, 0:ncol-1]
    image2 = image2[0:nrow-1, 0:ncol-1]
    
    r_scl=np.zeros(nrow)
    r_off=np.zeros(nrow)
    r_std=np.zeros(nrow)
    
    for jj in range(0,nrow-1):
        x=np.concatenate((np.ones((ncol-1,1)), np.expand_dims(image1[jj,].conj().transpose(), axis=1)), axis=1)#-1 to adjust to python indexing?
        y=image2[jj,].conj().transpose()
        bb, ab, aa, cc = np.linalg.lstsq(x,y)
        
        ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
        #ft=np.multiply(x[:,1]*bb[1]) + bb[0]

        adf=np.abs(np.squeeze(y)-np.squeeze(ft))
        sigma=np.std(np.squeeze(y)-np.squeeze(ft))

        inside = np.where(adf < 2*sigma)
        bb, ab, aa, cc = np.linalg.lstsq(x[inside[1],], y[inside[1]])
        
        ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
        
        r_scl[jj]=bb[1]
        r_off[jj]=bb[0]
        r_std[jj]=np.std(y[0]-ft[0])
        
    c_scl=np.zeros(nrow)
    c_off=np.zeros(nrow)
    c_std=np.zeros(nrow)
   
    for jj in range(0,ncol-1):
        
        x=np.concatenate((np.ones((nrow-1,1)), np.expand_dims(image1[:,jj], axis=1)), axis=1)
        y=image2[:,jj]
        bb, ab, aa, cc = np.linalg.lstsq(x,y)
        
        ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
        
        adf=np.abs(np.squeeze(y)-np.squeeze(ft))
        sigma=np.std(np.squeeze(y)-np.squeeze(ft))

        inside = np.where(adf < 2*sigma)
        bb, ab, aa, cc = np.linalg.lstsq(x[inside[1],], y[inside[1]])
            
        ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
        
        c_scl[jj]=bb[1]
        c_off[jj]=bb[0]
        c_std[jj]=np.std(y[0]-ft[0])
    
    nsz=(nrow-1)*(ncol-1)
    la_1=np.reshape(image1, (nsz,1))
    la_2=np.reshape(image2, (nsz,1))
    
    x=np.concatenate((np.ones((nsz,1)),la_1), axis=1)
    y=la_2
    bb, ab, aa, cc = np.linalg.lstsq(x,y)
    
    ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
    
    adf=np.abs(np.squeeze(y)-np.squeeze(ft))
    sigma=np.std(np.squeeze(y)-np.squeeze(ft))
    
    inside = np.where(adf < 2*sigma)
    bb, ab, aa, cc = np.linalg.lstsq(x[inside[1],], y[inside[1]])
    
        
    ft=np.squeeze([a*bb[1] for a in x[:,1]]) + bb[0]
    
    t_off=bb[0]
    t_scl=bb[1]
    t_std=np.std(y[0]-ft[0])
    
    rows=0
    
    return t_off, t_scl, t_std

def compensate_bad_columns(image, header):
    #LM 200127 This does not need to be used since it is already done in the OBC says Georgi.
    
    #this is a function to compensate bad columns if in the image
    
    ncol=int(header['NCOL'])+1
    nrow=int(header['NROW'])
    
    ncolskip=int(header['NCSKIP'])
    
    ncolbinC=int(header['NCBIN'])
    ncolbinF=2**int(header['NColBinFPGA'])
    
    #change to Leading if Trailing does not work properly
    blank=int(header['TBLNK'])
    
    gain=2**(int(header['DigGain']))
    
    if ncolbinC == 0: #no binning means binning of one
        ncolbinC=1
    
    if ncolbinF==0: #no binning means binning of one
        ncolbinF=1
    
    #bad column analysis
    
    n_read, n_coadd = binning_bc(ncol, ncolskip, ncolbinF, ncolbinC, np.asarray(header['BC']))
    
    if ncolbinC>1:
        for j_c in range(0,ncol):
            if ncolbinC*ncolbinF != n_coadd[j_c]:
                #remove gain adjustment
                image[0:nrow-1,j_c] = image[0:nrow-1,j_c]*gain
                
                #remove added superpixel value due to bad columns and read out offset
                image[0:nrow-1,j_c] = image[0:nrow-1,j_c] - n_read[j_c]*(blank-128) -128
                
                #multiply by number of binned column to actual number readout ratio
                image[0:nrow-1,j_c] = image[0:nrow-1,j_c] * ((ncolbinC*ncolbinF)/n_coadd[j_c])
                
                #add number of FPGA binned
                image[0:nrow-1,j_c] = image[0:nrow-1,j_c] + ncolbinF*(blank-128) + 128
                
                #add gain adjustment back
                image[0:nrow-1,j_c] = image[0:nrow-1,j_c]/gain
                
                print('Col: ',j_c,', n_read: ',n_read[j_c],', n_coadd: ',n_coadd[j_c],', binned pixels: ',ncolbinC*ncolbinF)
    
    return image

#
#def get_true_image_from_compensated(image, header):
#    
#    #calculate true image by removing readout offset, pixel blank value and
#    #normalising the signal level according to readout time
#    
#    #remove gain
#    true_image = image * 2**(int(header['DigGain']) & 255)
#    
#    for j_c in range(0,int(header['NCol'])):
#        true_image[0:header['NRow'], j_c] = ( true_image[0:header['NRow'],j_c] - 
#                  2**header['NColBinFPGA'] * (header['BlankTrailingValue']-128) - 128 )
#    
#    
#    return true_image





    
