#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 13:33:31 2020

@author: lindamegner
"""




from mats_l1_processing.items_units_functions import read_files_in_protocol_as_ItemsUnits
from mats_l1_processing.experimental_utils import  plot_CCDimage

from mats_l1_processing.experimental_utils import readprotocol
import matplotlib.pyplot as plt
import numpy as np

from L1_calibration_functions import  L1_calibrate_only_partly




#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/Diffusor/DiffusorFlatTests/'
#protocol='ABOUT.txt'

directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/LimbFlatfield/'
#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200330_flatfields_0C/'
protocol='protocol_special.txt'



#directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/Diffusor/200910_flatfields/'
#protocol='protocol_SensitivityToDiffusorGeometry.txt'

read_from='imgview'  
df_protocol=readprotocol(directory+protocol)

# df_bright=df_protocol[df_protocol.DarkBright=='B']
# CCDitems=read_all_files_in_protocol(df_bright, read_from,directory)

ItemsUnits=read_files_in_protocol_as_ItemsUnits(df_protocol,directory,3, read_from)



fig,ax=plt.subplots(6,1,figsize=(6,7.5))
fig2,ax2=plt.subplots(2,1,figsize=(8,5.5))
fig3,ax3=plt.subplots(2,1,figsize=(8.5,5))

ydiffdict={'IR1': 47,'IR2': 76,'IR3': 37,'IR4': 0,'UV1': 15,'UV2': 192 }#Values from Alignment.xls or AlignmentFromFitting.xls
xdiffdict={'IR1': -75,'IR2': 144,'IR3': 66,'IR4': 0,'UV1': 88,'UV2': 156}#Values from Alignment.xls or AlignmentFromFitting.xls


for ind,ItemsUnit in enumerate(ItemsUnits):

    channel=ItemsUnit.imageItem['channel']
    print(channel)

    img1=ItemsUnit.subpic
    #img1= L1_calibrate_only_partly(ItemsUnit.imageItem)
    
    
    
    mean_vertical=np.mean(img1[:,100:1948],1)
    mean_horizontal=np.mean(img1[300:400,:],0)   

    
    plot_CCDimage(img1, fig, ax[ind],title=channel)
    
    yaxis=np.arange(img1.shape[0])
    xaxis=np.arange(img1.shape[1])
    
    yaxisshifted=yaxis-ydiffdict[channel] 
    xaxisshifted=xaxis-xdiffdict[channel]
    
    ax2[0].plot(mean_vertical,yaxis, label=channel) 
    if ydiffdict[channel]>=-900:
        ax2[1].plot(mean_vertical,yaxisshifted, label=channel) 
        
    ax3[0].plot(xaxis, mean_horizontal,label=channel) 
    if xdiffdict[channel]>=-900:
        ax3[1].plot(xaxisshifted, mean_horizontal, label=channel)       

ax2[0].legend()
ax2[1].legend()
ax3[0].legend()
ax3[1].legend()



fig.suptitle('Flatfield Shape')
fig.savefig('Limb_flatfield_shape'+'.jpg')
fig2.suptitle('Mean vertical flatfield signal')
fig2.savefig('Limb_flatfield_vertical'+'.jpg')
fig3.suptitle('Mean horizontal flatfield signal')
fig3.savefig('Limb_flatfield_horizontal'+'.jpg')