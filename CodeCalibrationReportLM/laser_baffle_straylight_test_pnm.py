#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 10:27:10 2020

@author: lindamegner

Coompares two images from pnmfile 
"""

import os
import glob
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt

#from experimental_utils import plot_CCDimage


def plot_CCDimage(image, fig, axis, title="", clim=999):
    import numpy as np


    #sp = axis.pcolormesh(image, cmap=plt.cm.jet)
    sp = axis.imshow(image)# cmap=plt.cm.jet)
    if clim == 999:
        mean = np.mean(image)
        std = np.std(image)
        sp.set_clim([mean - 2 * std, mean + 2 * std])
    else:
        sp.set_clim(clim)
    fig.colorbar(sp, ax=axis)
    axis.set_title(title)
    return sp


def diffplot(img1,img2, fig,ax, title1, title2, suptitle):

    #img2_scaled=img2*img1.mean()/img2.mean()
    diff=img1-img2#_scaled
    

    plot_CCDimage(img1, fig, ax[0],title=title1)
    plot_CCDimage(img2, fig, ax[1],title=title2)
    plot_CCDimage(diff, fig, ax[2],title='difference panel1-panel2')
    
    fig.suptitle(suptitle)
    return fig
    
    #fig.savefig('Limb_baffle_flatfield_diff_'+when+'_'+ channel+'.jpg')


#directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210519StraylightTest/PayloadImages/'



directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210526OHBBaffleInFoVTest/PayloadImages/'


channellist=['IR1', 'IR2', 'IR3', 'IR4']
for channel in channellist:

        
    #Sequence 1: #Laser shining at last vane (inner vane) . Note that the light scatters into the instrument to the side of #the mirror. The vane is hit on the lower side, ie where Jörgs mountain chain was situated. Ie. The side #closest to the baffle mirror.
    if channel=='IR1':
        seq1=directory+'1306064591647781376_1'+'.pnm'
    elif channel=='IR2':
        seq1=directory+'1306064776599288832_2'+'.pnm'
    elif channel=='IR3':
        seq1=directory+'1306064878954131968_3'+'.pnm'
    elif channel=='IR4':
        seq1=directory+'1306064962398284800_4'+'.pnm'

        #Sequence 2: #Dark pictures
    if channel=='IR1':
        seq2=directory+'1306065064952361984_1'+'.pnm'
    elif channel=='IR2':
        seq2=directory+'1306065137242980864_2'+'.pnm'
    elif channel=='IR3':
        seq2=directory+'1306065211470046976_3'+'.pnm'
    elif channel=='IR4':
        seq2=directory+'1306065285403350784_4'+'.pnm'

        #Sequence 3: # Laser shining at first vane (outer vane) . The vane is hit on the lower side, ie where Jörgs mountain chain was situated. Ie. The side closest to the baffle mirror
    if channel=='IR1':
        seq3=directory+'1306065560305511424_1'+'.pnm'
    elif channel=='IR2':
        seq3=directory+'1306065629451827968_2'+'.pnm'
    elif channel=='IR3':
        seq3=directory+'1306065698850433280_3'+'.pnm'
    elif channel=='IR4':
        seq3=directory+'1306065769960586496_4'+'.pnm'

        #Sequence 4: #Upper side of baffle. Laser shining at last vane (inner vane) . Note that the light scatters into the #instrument to the side of the mirror. The vane is hit on the upper side, ie not where Jörgs mountain #chain was situated. Ie. The side furtest the baffle mirror.
    if channel=='IR1':
        seq4=directory+'1306068103546707200_1'+'.pnm'
    elif channel=='IR2':
        seq4=directory+'1306068177691040000_2'+'.pnm'
    elif channel=='IR3':
        seq4=directory+'1306068250506317056_3'+'.pnm'
    elif channel=='IR4':
        seq4=directory+'1306068347324310272_4'+'.pnm'

#Sequence 5: #Dark pictures
    if channel=='IR1':
        seq5=directory+'1306068808982605056_1'+'.pnm'
    elif channel=='IR2':
        seq5=directory+'1306068883717270016_2'+'.pnm'
    elif channel=='IR3':
        seq5=directory+'1306068961775573760_3'+'.pnm'
    elif channel=='IR4':
        seq5=directory+'1306069023806991616_4'+'.pnm'
        
#Sequence 6: # Laser shining at first vane (outer vane) . The vane is hit on the upper side, ie not where Jörgs mountain chain was situated. Ie. The side furthest the baffle mirror.
    if channel=='IR1':
        seq6=directory+'1306068527946090752_1'+'.pnm'
    elif channel=='IR2':
        seq6=directory+'1306068596482040320_2'+'.pnm'
    elif channel=='IR3':
        seq6=directory+'1306068666673431296_3'+'.pnm'
    elif channel=='IR4':
        seq6=directory+'1306068729693557760_4'+'.pnm'

    imgseq1 = np.float64(Image.open(seq1))  # read image
    imgseq2 = np.float64(Image.open(seq2))
    imgseq3 = np.float64(Image.open(seq3))  # read image
    imgseq4 = np.float64(Image.open(seq4))
    imgseq5 = np.float64(Image.open(seq5))  # read image
    imgseq6 = np.float64(Image.open(seq6))  # read image
    
    
    fig,ax=plt.subplots(3,1)
    sp=diffplot(imgseq6,imgseq5, fig,ax, title1='Laser at outer upper vane', title2='No laser', suptitle=channel)
    fig.savefig('LaserBaffeInViewOuterUpperVane'+ channel+'.jpg')

