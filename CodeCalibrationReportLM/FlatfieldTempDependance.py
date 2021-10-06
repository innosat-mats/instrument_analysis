#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:25:11 2020

@author: lindamegner
"""


from mats_l1_processing.LindasCalibrationFunctions import read_files_in_protocol_as_ItemsUnits
from mats_l1_processing.read_in_functions import readprotocol
import matplotlib.pyplot as plt
import numpy as np


directory='/Users/lindamegner/MATS/retrieval/Calibration/AfterLightLeakage/Flatfields/20200511_temperature_dependence/'



read_from='rac'  

df_protocol=readprotocol(directory+'protocol.txt')
#df_only2 = df_protocol[(df_protocol.index-2) % 3 != 0]

CCDItemsUnits=read_files_in_protocol_as_ItemsUnits(df_protocol,directory,3,read_from)

fig_all=plt.figure()
ax_all=fig_all.gca()

channels=['IR1','IR2','IR3','IR4','UV1','UV2']
for channel in channels:
    CCDItemsUnitsSelect= list(filter(lambda x: ( x.imageItem['channel']==channel),CCDItemsUnits))
    
    
    if channel[0]=='I':
        maxplot=8
    elif channel[0]=='U':
        maxplot=7
    else:
        raise Exception('Unknown channel name')
    fig,ax=plt.subplots(maxplot,1)
    imgmean=[]
    imgmeanerr=[]
    temperature=[]
    for i, ItemsUnit in enumerate(CCDItemsUnitsSelect):
         # fig=plt.figure()
         # ax=fig.gca()    
         ItemsUnit.plot(fig,ax[i], title=ItemsUnit.imageItem['channel'])
         ax[i].text(100,200,'mean:'+str(np.mean(ItemsUnit.subpic)))
         ax[i].text(100,320,'T:'+str(ItemsUnit.imageItem['temperature']))
         subpixmid=ItemsUnit.subpic[100:400,100:1948]
         imgmean.append(np.mean(subpixmid))
         imgmeanerr.append(np.std(subpixmid)/np.sqrt(subpixmid.size))
         temperature.append(ItemsUnit.imageItem['temperature'])
    fig.suptitle('Flat fields as function of temperature')
   # fig.savefig('FlatfeildTempVariability_200617_'+channel+'.jpg')
    

    imgmeanmean=np.mean(imgmean)
    ax_all.errorbar(temperature, imgmean/imgmean[0],yerr=imgmeanerr/imgmean[0], label=channel)

ax_all.legend()

fig_all.suptitle('Mean of flat field as func of temp')
fig_all.savefig('FlatfieldTempDependence.jpg')
