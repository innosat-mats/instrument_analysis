#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 09:34:19 2019

@author: franzk
"""

import numpy as np
import matplotlib.pyplot as plt

# Use from L1_functions 
#from L1_functions import readimgpath, predict_image, get_true_image, desmear_true_image, compensate_bad_columns, get_true_image_from_compensated



# Use from L1_calibration_functions
from mats_l1_processing.L1_calibration_functions import  get_true_image, desmear_true_image, compensate_bad_columns,predict_image, readimgpath, get_true_image_from_compensated




def add_capital_naming_to_header(CCDitem):
    try:
        CCDitem["NCOL"]
    except:
        CCDitem["NCOL"] = CCDitem["NCol"] 
    try:
        CCDitem["NROW"]
    except:
        CCDitem["NROW"] = CCDitem["NRow"] 
    try:
        CCDitem["NCSKIP"]
    except:
        CCDitem["NCSKIP"] = CCDitem["NColSkip"] 

    try:
        CCDitem["TBLNK"]
    except:
        CCDitem["TBLNK"] = CCDitem["BlankTrailingValue"]
    try:
        CCDitem["LBLNK"]
    except:
        CCDitem["LBLNK"] = CCDitem["BlankLeadingValue"]
    try:
        CCDitem["NRBIN"]
    except:
        CCDitem["NRBIN"] = CCDitem["NRowBinCCD"] 
    try:
        CCDitem["NCBIN"]
    except:
        CCDitem["NCBIN"] = CCDitem["NCBIN CCDColumns"] #Changed fron NColBinCCD 220908 by LM
    try:
        CCDitem["NRSKIP"]
    except:
        CCDitem["NRSKIP"] = CCDitem["NRowSkip"] 
    try:
        CCDitem["NCSKIP"]
    except:
        CCDitem["NCSKIP"] = CCDitem["NColSkip"] 
        
    try:
        CCDitem["NFLUSH"]
    except:
        CCDitem["NFLUSH"] = CCDitem["N_flush"] 
                
    try:
        CCDitem["TEXPMS"]
    except:
        CCDitem["TEXPMS"] = CCDitem["Texposure"] 
              

    try:
        CCDitem["DigGain"] = CCDitem["GAIN Truncation"]
    except:
        CCDitem["DigGain"] = CCDitem["GAIN"] & 0b1111
    #Used to be the following but I think that is wrong LM220909: CCDitem["DigGain"] = CCDitem["Gain"]
    
    
    try:
        CCDitem["BC"]
    except:
        CCDitem["BC"] = CCDitem["BadCol"]
    #CCDitem["BC"] =np.array(int(CCDitem["BC"]))

    #Convert BC to list of integer instead of str      
    if CCDitem['BC']=='[]':
        CCDitem['BC'] =np.array([])
    else:
        # strlist=CCDitem['BC'][1:-1].split(' ')
        # CCDitem['BC'] = np.array([int(i) for i in strlist])
        
        CCDitem["BC"] =[int(x) for x in CCDitem["BC"]]

    
get_true_from_L1_calibration_functions=True

# MATS data chain for raw images
# this script deals with:
# 1. Image prediction from reference image according to header
# 2. Applying bad column correction and removing offsets to get "true image"
# 3. Compensating bad columns in MATS payload OBC


# Import two reference images for the CCD in high and low signal modes
# Reference image depends on the temperature


ref_hsm_image, hsm_header = readimgpath('/Users/lindamegner/MATS/retrieval/Level1/data/2019-02-08 rand6/', 0, 0)
ref_lsm_image, lsm_header = readimgpath('/Users/lindamegner/MATS/retrieval/Level1/data/2019-02-08 rand6/', 4, 0)

recv_image, recv_header = readimgpath('/Users/lindamegner/MATS/retrieval/Level1/data/2019-02-08 rand6/', 32, 0) # Image nr 32 was used to check bc too

#recv_image, recv_header = readimgpath('/Users/lindamegner/MATS/retrieval/Level1/data/2019-02-08 rand6/', 30, 0)

# =============================================================================
# Name of directory on Box: 2019-02-08_rand_Georgi.xlsx
# =============================================================================


add_capital_naming_to_header(hsm_header)
add_capital_naming_to_header(lsm_header)
add_capital_naming_to_header(recv_header)


#Step 1
# Predict the received image from reference image according to header
pred_image, pred_header = predict_image(ref_hsm_image.copy(), hsm_header.copy(), ref_lsm_image.copy(), lsm_header.copy(), recv_header.copy())

add_capital_naming_to_header(pred_header)
    

    
#Step 2
# Make actual image out both received and predicted images
# These can be compared if no compression is used.
    


recv_true_image = get_true_image(recv_header.copy(), recv_image.copy())

recv_true_image = desmear_true_image(recv_header.copy(), recv_true_image.copy())

pred_true_image = get_true_image(pred_header.copy(), pred_image.copy())

#Step 3
# Bad column compensation for MATS payload OBC
recv_comp_image = compensate_bad_columns(recv_header.copy(), recv_image.copy())


    
#Step 4
#Getting "true image" from compensated one in MATS payload OBC
    
true_comp_image = get_true_image_from_compensated(recv_comp_image.copy(), recv_header.copy())


if get_true_from_L1_calibration_functions: 
    true_comp_image = desmear_true_image(recv_header.copy(), true_comp_image.copy())
else: # get_true_from_L1_functions
    true_comp_image = desmear_true_image(true_comp_image.copy(), recv_header.copy())

print('recv_true_img',recv_true_image[-4:,-4:])
# =============================================================================
#Lindas plotting
# mean_img = np.mean(np.mean(recv_image))
# clim=[mean_img-image_display_adjustment, mean_img+image_display_adjustment]
# fig1=diffplot(recv_image, pred_image,'CCD image','predicted image',clim=clim)
# 
# 
# true_mean_img = np.mean(np.mean(recv_true_image))
# clim=[true_mean_img-image_display_adjustment, true_mean_img+image_display_adjustment]
# fig2=diffplot(recv_true_image, pred_true_image,'CCD true image','predicted true image',clim=clim)
# 
# fig3=diffplot(recv_comp_image, true_comp_image,'compensated in OBC image','True image after compensation in OBC software',clim=clim)
# 
# =============================================================================


#  plotting
mean_img = np.mean(np.mean(recv_image))
std_img = np.std(recv_image)

image_display_adjustment = 600#2*std_img#600
image_display_adjustment2 = 100




def georgiplot(fig, ax,image,vmin,vmax, title):
    sp=ax.imshow(image, vmin=vmin, vmax=vmax, aspect='auto',origin='lower', interpolation='none',cmap='jet')   
    #xaxis=range(0, image.shape[1])
    #yaxis=range(0,image.shape[0])
    #sp=ax.pcolormesh(xaxis, yaxis, image, vmin=vmin, vmax=vmax, shading='neaest',cmap='jet')
    
    plt.title(title)
    ax.set_xlabel('Pixels')
    ax.set_ylabel('Pixels')
    ax.set_title(title)
    fig.colorbar(sp, ax=ax)
    return sp


fig1, ax=plt.subplots(3,1)
sp=georgiplot(fig1, ax[0], recv_image, vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment, title='CCD image')
sp=georgiplot(fig1, ax[1], pred_image, vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment, title='predicted image')
mean_diff=np.mean(np.mean(recv_image-pred_image))
sp=georgiplot(fig1, ax[2],recv_image-pred_image, vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, title='diff image')
fig1.savefig('GeorgisModelPredictedVersusReceived.jpg')

true_mean_img = np.mean(np.mean(recv_true_image))
fig1, ax=plt.subplots(3,1)
sp=georgiplot(fig1, ax[0],recv_true_image, vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment, title= 'CCD true image')
sp=georgiplot(fig1, ax[1],pred_true_image, vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment, title= 'predicted true image')
mean_diff=np.mean(np.mean(recv_true_image-pred_true_image))
sp=georgiplot(fig1, ax[2],recv_true_image-pred_true_image,  vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, title='diff image' )
fig1.savefig('GeorgisModelPredictedTrueVersusReceivedTrue.jpg')

fig1, ax=plt.subplots(3,1)
sp=georgiplot(fig1, ax[0],recv_comp_image, vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment, title='compensated in OBC image')
sp=georgiplot(fig1, ax[1],true_comp_image, vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment, title='True image after compensation in OBC software')
mean_diff=np.mean(np.mean(recv_comp_image-true_comp_image))
sp=georgiplot(fig1, ax[2],recv_comp_image-true_comp_image, vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, title='diff image')
fig1.savefig('GeorgisModelTrueCompensatedVersusReceivedCompensated.jpg')


""" Georgis origingal code below
fig1, ax=plt.subplots(3,1)
plt.subplot(ax[0])
plt.imshow(recv_image, cmap='jet', vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment,  aspect='auto')
plt.title('CCD image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()


plt.subplot(ax[1])
plt.imshow(pred_image, cmap='jet', vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment, aspect='auto')
plt.title('predicted image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

mean_diff=np.mean(np.mean(recv_image-pred_image))
plt.subplot(ax[2])
plt.imshow(recv_image-pred_image, cmap='jet', vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, aspect='auto')
#plt.imshow(recv_image-pred_image, cmap='jet', vmin=-50, vmax=300,  extent=[0,75,0,511], aspect='auto')

plt.title('diff image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

fig1.savefig('GeorgisModelPredictedVersusReceived.jpg')

true_mean_img = np.mean(np.mean(recv_true_image))

fig1, ax=plt.subplots(3,1)
plt.subplot(ax[0])
plt.imshow(recv_true_image, cmap='jet', vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment, aspect='auto')
plt.title('CCD true image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

plt.subplot(ax[1])
plt.imshow(pred_true_image, cmap='jet', vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment, aspect='auto')
plt.title('predicted true image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

mean_diff=np.mean(np.mean(recv_true_image-pred_true_image))
plt.subplot(ax[2])
plt.imshow(recv_true_image-pred_true_image, cmap='jet', vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, aspect='auto')
plt.title('diff image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

fig1.savefig('GeorgisModelPredictedTrueVersusReceivedTrue.jpg')

fig1, ax=plt.subplots(3,1)
plt.subplot(ax[0])
plt.imshow(recv_comp_image, cmap='jet', vmin=mean_img-image_display_adjustment, vmax=mean_img+image_display_adjustment, aspect='auto')
plt.title('compensated in OBC image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

plt.subplot(ax[1])
plt.imshow(true_comp_image, cmap='jet', vmin=true_mean_img-image_display_adjustment, vmax=true_mean_img+image_display_adjustment,  aspect='auto')
plt.title('True image after compensation in OBC software')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()

mean_diff=np.mean(np.mean(recv_comp_image-true_comp_image))
plt.subplot(ax[2])
plt.imshow(recv_comp_image-true_comp_image, cmap='jet', vmin=mean_diff-image_display_adjustment2, vmax=mean_diff+image_display_adjustment2, aspect='auto')
plt.title('diff image')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.colorbar()
plt.show()


fig1.savefig('GeorgisModelTrueCompensatedVersusReceivedCompensated.jpg')
#, extent=[0,75,0,511]

"""