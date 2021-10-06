#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 07:28:46 2020

Forward model for MATS' calibration

@author: lindamegner
"""

import numpy as np
import matplotlib.pyplot as plt

from mats_l1_processing.LindasCalibrationFunctions import read_all_files_in_protocol
from mats_l1_processing.read_in_functions import readprotocol
from mats_l1_processing.L1_calibration_functions import get_true_image, desmear_true_image, subtract_dark, compensate_flatfield, CCD
from mats_l1_processing.L1_calibration_functions import calculate_flatfield, calculate_dark, desmear_true_image_reverse, get_true_image_reverse, bin_image_using_predict_and_get_true_image, bin_image_with_BC 

#from LindasCalibrationFunctions import plot_CCDimage 




def plot(pic,fig,axis,title='',clim=999):

                    
    sp=axis.pcolormesh(pic,cmap=plt.cm.jet)
    axis.set_title(title)
    if clim==999:
        mean=pic.mean()
        std=pic.std()
        sp.set_clim([mean-1*std,mean+1*std])
       
    else:
        sp.set_clim(clim)

    fig.colorbar(sp,ax=axis)

    return sp    
    


def forward(photons,CCDitem):
    simage_raw=np.float64(photons*np.ones([511,2048]))
    #simage_raw=np.float64(photons*np.ones_like(CCDitem['IMAGE']))
    # TOD0 Step 8 Transform from photons to electrons and then to LSB.   

    #if CCDitem['NCBIN CCDColumns']>1 or CCDitem['NCBIN FPGAColumns']>1 or CCDitem['NRBIN']>1 : 
    simage_raw_binned=bin_image_using_predict_and_get_true_image(CCDitem.copy(), simage_raw.copy())
   # simage_raw_binned=get_true_image(CCDitem.copy(),simage_raw_binned)  
    
    #plotmean=photons*CCDitem['NCBIN CCDColumns']*CCDitem['NCBIN FPGAColumns']*CCDitem['NRBIN']
    #clims=[plotmean-np.sqrt(plotmean), plotmean+np.sqrt(plotmean)]
    plot(simage_raw_binned, fig, ax[0,f], title='raw simulated image')
    
    # Now modify the image in forward direction and Plot the result
    
    
    # Step 7 Add ghost imaging. TBD.
    
    
    #flatfield
    CCDitem_nobin=CCDitem.copy()
    CCDitem_nobin['NCBIN CCDColumns']=1 
    CCDitem_nobin['NCBIN FPGAColumns']=1
    CCDitem_nobin['NRBIN']=1    
    image_flatf_fact=calculate_flatfield(CCDitem_nobin.copy())      
    simage_flatf=simage_raw*image_flatf_fact 
    #simage_flatf=simage_raw*image_flatf_fact[CCDitem['NRSKIP']:CCDitem['NRSKIP']+CCDitem['NROW'],
    #                                     CCDitem['NCSKIP']:CCDitem['NCSKIP']+CCDitem['NCOL']+1]


    simage_flatf_binned=bin_image_using_predict_and_get_true_image(CCDitem.copy(), simage_flatf.copy())
    #simage_flatf_binned=get_true_image(CCDitem.copy(),simage_flatf_binned)
    plot(simage_flatf_binned, fig, ax[1,f], title='raw+flat')
    
    #dark
    # TBD: Decide on threshold fro when to use pixel correction (little dark current) and when to use average image correction (large dark current). 
    # TBD: The temperature needs to be decided in a better way then taken from the ADC as below.
    # Either read from rac files of temperature sensors or estimated from the top of the image
    
    simage_dark=simage_flatf+calculate_dark(CCDitem_nobin.copy())
    #dark_fullpic=calculate_dark(CCDitem.copy())    
    #simage_dark=simage_flatf+dark_fullpic[CCDitem['NRSKIP']:CCDitem['NRSKIP']+CCDitem['NROW'],
    #                                      CCDitem['NCSKIP']:CCDitem['NCSKIP']+CCDitem['NCOL']+1]

    
    
    
    
    # Binning and bad columns
    simage_dark_binned=bin_image_using_predict_and_get_true_image(CCDitem.copy(), simage_dark.copy())
    # testfig=plt.figure()
    # testax=testfig.gca()
    # hej=simage_binned-simage_dark
    # plot(hej, testfig, testax, title='in-out')
   
    #simage_dark_binned=get_true_image(CCDitem.copy(),simage_dark_binned)
    
    
    plot(simage_dark_binned, fig, ax[2,f], title='raw+flat+dark+binned')
    
    #add smear
    simage_smear=desmear_true_image_reverse(CCDitem.copy(), simage_dark_binned.copy())
    plot(simage_smear, fig, ax[3,f], title='raw+flat+dark+smear')
    
    
        
    #add bias
    simage_bias=get_true_image_reverse(CCDitem.copy(),simage_smear.copy())
    plot(simage_bias,fig, ax[4,f], title='raw+flat+dark+smear+bias')
    
    return simage_raw,simage_raw_binned, simage_flatf, simage_flatf_binned, simage_dark, simage_dark_binned, simage_smear, simage_bias
        
    
    
    

def backward(input_image,CCDitem):    
    
        # # Do normal calibration to reverse the forward model. Same as in L1_calibrate but plots between the different steps. 
    image=input_image.copy()
    plot(image,fig, ax[4,b], 'From forward')
    plot(input_image-image,fig, ax[4,d], 'simage-image')

    image_bias_sub = get_true_image(CCDitem, image)
    plot(image_bias_sub,fig, ax[3,b], 'Bias subtracted') 
    plot(simage_smear-image_bias_sub,fig, ax[3,d], 'simage-image')
    
    image_desmeared = desmear_true_image(CCDitem,image_bias_sub.copy())
    plot(image_desmeared,fig, ax[2,b],' Desmeared LSB')  
    plot(simage_dark_binned-image_desmeared,fig, ax[2,d], 'simage-image')
    
    image_dark_sub=subtract_dark(CCDitem,image_desmeared.copy())
    plot(image_dark_sub,fig, ax[1,b], ' Dark current subtracted.')     
    plot(simage_flatf_binned-image_dark_sub,fig, ax[1,d], 'simage-image')
    
    
    image_flatf_comp=compensate_flatfield(CCDitem,image_dark_sub.copy())
    #plotmean=photons*CCDitem['NCBIN CCDColumns']*CCDitem['NCBIN FPGAColumns']*CCDitem['NRBIN']
    #clims=[plotmean-np.sqrt(plotmean), plotmean+np.sqrt(plotmean)]

    plot(image_flatf_comp,fig, ax[0,b], ' Flat field compensated.')     
    plot(simage_raw_binned-image_flatf_comp,fig, ax[0,d], 'simage-image')
    
    fig.suptitle('Forward model followed by backward i.e calibration')
    

    return image, image_bias_sub, image_desmeared, image_dark_sub, image_flatf_comp



# =============================================================================
# Main
# =============================================================================




clims=[-2,2]

# Read in a CCDitem 



#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/Diffusor/DiffusorFlatTests/'
#protocol='ForwardModelTestProto.txt'


directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Binning/Binning-simulation/'
protocol='PROTOCOL-BINNING.txt'

read_from='rac'  
df_protocol=readprotocol(directory+protocol)
df_bright=df_protocol[df_protocol.DarkBright=='B']
CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)
# The imagge of this CCDitem is not used  , only the meta data


CCDunits={}


CCDitem=CCDitems[0]

calibration_toml_file='/Users/lindamegner/MATS/retrieval/git/MATS-L1-processing/scripts/calibration_data_linda.toml'

#Check  if the CCDunit has been created. It takes time to create it so it should not be created if not needed
try: CCDitem['CCDunit']
except: 
    try:
        CCDunits[CCDitem['channel']]
    except:  
        CCDunits[CCDitem['channel']]=CCD(CCDitem['channel'], calibration_toml_file) 
    CCDitem['CCDunit']=CCDunits[CCDitem['channel']]


    

# OleMartin:  The below should happen in Level 0. /LM
CCDitem['BC']=np.array(CCDitem['BC'])  
    



fig,ax=plt.subplots(5,3)
f=0
b=1
d=2

# Create simage
photons= 1000

simage_raw,simage_raw_binned, simage_flatf, simage_flatf_binned, simage_dark, simage_dark_binned, simage_smear, simage_bias=forward(photons,CCDitem)

image, image_bias_sub, image_desmeared, image_dark_sub, image_flatf_comp=backward(simage_bias,CCDitem)

fig.savefig('forward_model.jpg')