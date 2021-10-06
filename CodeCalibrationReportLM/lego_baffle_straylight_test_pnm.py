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

#from LindasCalibrationFunctions import plot_CCDimage


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


directory='/Users/lindamegner/MATS/retrieval/Calibration/Final_AIT_2021/20210519StraylightTest/PayloadImages/'






channellist=['IR1', 'IR2', 'IR3', 'IR4', 'UV1', 'UV2']
for channel in channellist:

        
    #Sequence 1: 
    if channel=='IR1':
        seq1=directory+'1305449083882019072_1'+'.pnm'
    elif channel=='IR2':
        seq1=directory+'1305449192406860288_2'+'.pnm'
    elif channel=='IR3':
        seq1=directory+'1305449253595596288_3'+'.pnm'
    elif channel=='IR4':
        seq1=directory+'1305449304686126592_4'+'.pnm'
    elif channel=='UV1':
        seq1=directory+'1305449386754669312_5'+'.pnm'
    elif channel=='UV2':
        seq1=directory+'1305449524338348288_6'+'.pnm'
        #Sequence 2: 
    if channel=='IR1':
        seq2=directory+'1305450067208236800_1'+'.pnm'
    elif channel=='IR2':
        seq2=directory+'1305450115931365888_2'+'.pnm'
    elif channel=='IR3':
        seq2=directory+'1305450155696319488_3'+'.pnm'
    elif channel=='IR4':
        seq2=directory+'1305450194705871616_4'+'.pnm'
    elif channel=='UV1':
        seq2=directory+'1305450242514709504_5'+'.pnm'
    elif channel=='UV2':
        seq2=directory+'1305450332483245824_6'+'.pnm'
        #Sequence 3: 
    if channel=='IR1':
        seq3=directory+'1305450451349762048_1'+'.pnm'
    elif channel=='IR2':
        seq3=directory+'1305450497036560128_2'+'.pnm'
    elif channel=='IR3':
        seq3=directory+'1305450538885284352_3'+'.pnm'
    elif channel=='IR4':
        seq3=directory+'1305450590093139712_4'+'.pnm'
    elif channel=='UV1':
        seq3=directory+'1305450632752304128_5'+'.pnm'
    elif channel=='UV2':
        seq3=directory+'1305450718275024384_6'+'.pnm'
        #Sequence 4: 
    if channel=='IR1':
        seq4=directory+'1305450930883346432_1'+'.pnm'
    elif channel=='IR2':
        seq4=directory+'1305451110623245312_2'+'.pnm'
    elif channel=='IR3':
        seq4=directory+'1305451171598846464_3'+'.pnm'
    elif channel=='IR4':
        seq4=directory+'1305451225164352512_4'+'.pnm'
    elif channel=='UV1':
        seq4=directory+'1305451277812850944_5'+'.pnm'
    elif channel=='UV2':
        seq4=directory+'1305451364743850752_6'+'.pnm'
#Sequence 5: 
    if channel=='IR1':
        seq5=directory+'1305451731109222400_1'+'.pnm'
    elif channel=='IR2':
        seq5=directory+'1305451793113601792_2'+'.pnm'
    elif channel=='IR3':
        seq5=directory+'1305451842919067392_3'+'.pnm'
    elif channel=='IR4':
        seq5=directory+'1305451884125961216_4'+'.pnm'
    elif channel=='UV1':
        seq5=directory+'1305451930034545920_5'+'.pnm'
    elif channel=='UV2':
        seq5=directory+'1305452026729202176_6'+'.pnm'


    imgseq1 = np.float64(Image.open(seq1))  # read image
    imgseq2 = np.float64(Image.open(seq2))
    imgseq3 = np.float64(Image.open(seq3))  # read image
    imgseq4 = np.float64(Image.open(seq4))
    imgseq5 = np.float64(Image.open(seq5))  # read image

    
    
    fig,ax=plt.subplots(3,1)
    sp=diffplot(imgseq5,imgseq3, fig,ax, title1='Lego rod at edge of baffle futher in (seq 5)', title2='No lego rod (seq 3)', suptitle=channel)
    fig.savefig('Legotest_seq_5vs3'+ channel+'.jpg')

