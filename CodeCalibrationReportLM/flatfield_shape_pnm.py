#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:27:10 2020

@author: lindamegner

Coompares two images from pnmfile 
"""
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt


from mats_l1_processing.LindasCalibrationFunctions import plot_CCDimage

# def plot_CCDimage(image, fig, axis, title="", clim=999):
#     import numpy as np
#     import matplotlib.pyplot as plt

#     #sp = axis.pcolormesh(image, cmap=plt.cm.jet)
#     sp = axis.imshow(image, cmap=plt.cm.jet)
#     if clim == 999:
#         mean = np.mean(image)
#         std = np.std(image)
#         sp.set_clim([mean - 1 * std, mean + 1 * std])
#     else:
#         sp.set_clim(clim)
        
        
#     fig.colorbar(sp, ax=axis)
#     axis.set_title(title)
#     return sp


# def imshowimage(fig, ax, image):
#     plt.imshow(pic)
#     if clim == 999:
#         mean = pic.mean()
#         std = pic.std()
#         plt.clim([mean - 2 * std, mean + 2 * std])
#     else:
#         plt.clim(clim)
#     plt.colorbar()

channellist=['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']


fig,ax=plt.subplots(6,1,figsize=(6,7.5))
fig2,ax2=plt.subplots(2,1,figsize=(8,5.5))
fig3,ax3=plt.subplots(2,1,figsize=(8.5,5))

ydiffdict={'IR1': 47,'IR2': 76,'IR3': -999,'IR4': 0,'UV1': 15,'UV2': 192 }#Values from Alignment.xls or AlignmentFromFitting.xls
xdiffdict={'IR1': -75,'IR2': 144,'IR3': -999,'IR4': 0,'UV1': 88,'UV2': 156}#Values from Alignment.xls or AlignmentFromFitting.xls

for ind, channel in enumerate(channellist):

        
    #dir1='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/PayloadImages/'
    dir1='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/PayloadImages20210421_16_18/'
# =============================================================================
#     Note! In box the folder is named Flatfield20210421
# =============================================================================
        
    
    
    if channel=='IR1':
        filename1=dir1+'1303052550753753600_1'+'.pnm'
        dfilename1=dir1+'1303052736716888320_1'+'.pnm'
    elif channel=='IR2':
        filename1=dir1+'1303053042442397952_4'+'.pnm'
        dfilename1=dir1+'1303053108130874624_4'+'.pnm'
    elif channel=='IR3':
        filename1=dir1+'1303053614869430528_3'+'.pnm'
        dfilename1=dir1+'1303053680075119104_3'+'.pnm'
    elif channel=='IR4':
        filename1=dir1+'1303053921036712704_2'+'.pnm'
        dfilename1=dir1+'1303053986963058432_2'+'.pnm'
    elif channel=='UV2':
        filename1=dir1+'1303054301145446656_6'+'.pnm'
        dfilename1=dir1+'1303054436714050304_6'+'.pnm'
    elif channel=='UV1':
        filename1=dir1+'1303055407319107072_5'+'.pnm'
        dfilename1=dir1+'1303055745598846464_5'+'.pnm'    
    
    


    pic1 = np.float64(Image.open(filename1))  # read image
    picd1 = np.float64(Image.open(dfilename1))  # read image

    img1=pic1-picd1
    mean_vertical=np.mean(img1[:,100:2000],1)
    mean_horizontal=np.mean(img1[400:450,:],0)   

    
    plot_CCDimage(img1, fig, ax[ind],title=channel)
    
    yaxis=np.arange(img1.shape[0])
    xaxis=np.arange(img1.shape[1])
    
    yaxisshifted=yaxis-ydiffdict[channel] 
    xaxisshifted=xaxis-xdiffdict[channel]
    
    ax2[0].plot(mean_vertical,-yaxis, label=channel) 
    if ydiffdict[channel]>=-900:
        ax2[1].plot(mean_vertical,-yaxisshifted, label=channel) 
        
    ax3[0].plot(xaxis, mean_horizontal,label=channel) 
    if xdiffdict[channel]>=-900:
        ax3[1].plot(xaxisshifted, mean_horizontal, label=channel)         

ax2[0].legend()
ax2[1].legend()  
ax3[0].legend()
ax3[1].legend()
fig.suptitle('Flatfield Shape')
fig.savefig('Limb_flatfield_shape'+'.jpg')
fig2.savefig('Limb_flatfield_vertical'+'.jpg')
fig3.savefig('Limb_flatfield_horizontal'+'.jpg')